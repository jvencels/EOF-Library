#!/bin/bash

dockerImageName="eoflibrary/eof_elmer84_of6"

echo "$DOCKER_PASSWORD" | docker login -u "$DOCKER_USERNAME" --password-stdin

if [ "$TRAVIS_BRANCH" == "master" ]; then
  docker tag test_image $dockerImageName:latest
  docker push $dockerImageName:latest
fi

if [ -n "${TRAVIS_TAG}" ] ; then
  echo "Git commit tag name: ${TRAVIS_TAG}"
  docker tag test_image $dockerImageName:$TRAVIS_TAG
  docker push $dockerImageName:$TRAVIS_TAG

  docker tag test_image $dockerImageName:latest
  docker push $dockerImageName:latest
fi
