# decentralized-row-means
This repository contains code for decentralized row means written for the new COINSTAC simulator (`v4.0.0`).

## Files
1. local.py - computes the local row sums and number of columns from local data and sends results to remote site.
2. remote.py - aggreagates the local row sums and number of columns sent by local sites and returns the global row means.
3. compspec.json - computation specifications

## Written For
- Python 3.6.6
- coinstac-simulator 4.0.0

## Usage
1. Update `npm`:\
`sudo npm i npm@latest -g`
2. Install `coinstac-simulator`:\
`sudo npm i -g coinstac-simulator@4.0.0`
3. Clone this repository:\
`git clone https://github.com/rsilva8/decentralized-row-means`
4. Change directory:\
`cd decentralized-row-means`
5. Build the docker image (Docker must be running):\
`docker build -t decentralized-row-means .`
7. Run the code:\
`sudo coinstac-simulator`
