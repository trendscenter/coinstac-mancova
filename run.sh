sudo rm -r test/output -v;
echo "Building docker file...";
sudo docker build  -t coinstacteam/dmancova .;
echo "Running Simulator...";
sudo coinstac-simulator;
