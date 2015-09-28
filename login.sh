#!/bin/sh

curl -H "Content-Type: application/json" \
  -X POST \
  -k \
  -d  "{email: 'Administrator', password: 'Administrator'}" \
  https://localhost:8443/rest/accesstoken | \
  jq -r '.sessionId'
