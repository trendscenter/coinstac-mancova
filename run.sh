sudo rm -r test/output -v;
sudo rm -r test/transfer -v;
sudo rm -r test/input/remote/simulatorRun/* -v;
echo "Building docker file...";
sudo docker build  -t dmancova .;
echo "Running Simulator...";
sudo coinstac-simulator --silly;
