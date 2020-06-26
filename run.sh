sudo rm -r test/output -v;
echo "Building docker file...";
sudo docker build  -t dmancova .;
echo "Running Simulator...";
sudo coinstac-simulator --silly;
