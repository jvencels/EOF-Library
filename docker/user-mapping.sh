#!/bin/bash

# docker: mapping host uid and gid to user inisde container
# --------------------------------------------------
# Original file obtained from:
# https://gist.github.com/renzok/29c9e5744f1dffa392cf

# if both not set we do not need to do anything
if [ -z "${HOST_USER_ID}" -o -z "${HOST_USER_GID}" -o -z "${USER}" ]; then
    echo "Flags '-e HOST_USER_ID=\$(id -u) -e HOST_USER_GID=\$(id -g)' were not set. You might have permission problems."
    /bin/bash
else
  # reset user_?id to either new id or if empty old (still one of above
  # might not be set)
  USER_ID=${HOST_USER_ID:=$USER_ID}
  USER_GID=${HOST_USER_GID:=$USER_GID}

  LINE=$(grep -F "${USER}" /etc/passwd)
  # replace all ':' with a space and create array
  array=( ${LINE//:/ } )

  # home is 5th element
  USER_HOME=${array[4]}

  sed -i -e "s/^${USER}:\([^:]*\):[0-9]*:[0-9]*/${USER}:\1:${USER_ID}:${USER_GID}/"  /etc/passwd
  sed -i -e "s/^${USER}:\([^:]*\):[0-9]*/${USER}:\1:${USER_GID}/"  /etc/group

  chown -R ${USER_ID}:${USER_GID} ${USER_HOME}

  exec su - "${USER}"
fi
