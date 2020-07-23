FROM coinstacteam/coinstac-base:python3.7-buster-beta
ENV MCRROOT=/usr/local/MATLAB/MATLAB_Runtime/v91
ENV MCR_CACHE_ROOT=/tmp

RUN apt-get clean && apt-get update && apt-get install -y \
    zip unzip wget \
    libx11-dev libxcomposite-dev \
    libxcursor-dev libxdamage-dev libxext-dev \
    libxfixes-dev libxft-dev libxi-dev \
    libxrandr-dev libxt-dev libxtst-dev \
    libxxf86vm-dev libasound2-dev libatk1.0-dev \
    libcairo2-dev gconf2 \
    libsndfile1-dev libxcb1-dev libxslt-dev \
    curl \
    libgtk-3-dev 

RUN mkdir /tmp/mcr_installer && \
    cd /tmp/mcr_installer && \
    wget http://ssd.mathworks.com/supportfiles/downloads/R2016b/deployment_files/R2016b/installers/glnxa64/MCR_R2016b_glnxa64_installer.zip && \
    unzip MCR_R2016b_glnxa64_installer.zip && \
    ./install -mode silent -agreeToLicense yes && \
    rm -Rf /tmp/mcr_installer

# Copy the current directory contents into the container
#WORKDIR /app
COPY requirements.txt /app

# Install any needed packages specified in requirements.txt
#RUN pip install -r requirements.txt
COPY ./groupicatv4.0b/icatb/nipype-0.10.0/nipype/interfaces/gift /usr/local/lib/python3.6/site-packages/nipype/interfaces/gift
RUN chmod -R a+wrx /app
#RUN chmod -R a+wrx /usr/local/MATLAB/MATLAB_Runtime/v91

COPY . /app
#RUN (timeout 300 bash /app/groupicatv4.0b/GroupICATv4.0b_standalone/run_groupica.sh /usr/local/MATLAB/MATLAB_Runtime/v91; exit 0)

ENV MCRROOT=/usr/local/MATLAB/MATLAB_Runtime/v91
ENV MCR_CACHE_ROOT=/computation/mcrcache

# Copy the current directory contents into the container
WORKDIR /computation
COPY requirements.txt /computation

# Install any needed packages specified in requirements.txt
RUN pip install -r requirements.txt
RUN pip install awscli s3utils
RUN pip install nipy

RUN mkdir -p /computation/mcrcache

RUN mkdir /output

RUN (timeout 20s /app/groupicatv4.0b/GroupICATv4.0b_standalone/run_groupica.sh /usr/local/MATLAB/MATLAB_Runtime/v91/; exit 0)

COPY ./coinstac_masking /computation/coinstac_masking
COPY ./coinstac_decentralized_row_means /computation/coinstac_decentralized_row_means
COPY ./groupicatv4.0b /computation/groupicatv4.0b
COPY ./coinstac_node_ops /computation/coinstac_node_ops
COPY ./coinstac_spatially_constrained_ica /computation/coinstac_spatially_constrained_ica
COPY ./local_data /computation/local_data

COPY ./*.py /computation/


COPY ./coinstac_mancova /computation/coinstac_mancova
