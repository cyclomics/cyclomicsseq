
CONTAINER_VERSION="0.6.0"
CONTAINER_NAME="cycas"

docker build . -t cyclomics/$CONTAINER_NAME:$CONTAINER_VERSION
docker tag cyclomics/$CONTAINER_NAME:$CONTAINER_VERSION cyclomics/$CONTAINER_NAME:latest
docker push cyclomics/$CONTAINER_NAME:$CONTAINER_VERSION
docker push cyclomics/$CONTAINER_NAME:latest

