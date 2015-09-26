#!/bin/sh

OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
file=""
sessionid=""
entryid=""

while getopts "e:s:f:" opt; do
    case "$opt" in
    e)  entryid=$OPTARG
        ;;
    s)  sessionid=$OPTARG
        ;;
    f)  file=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

curl --form "file=@${file}" \
  -H "X-ICE-Authentication-SessionId: ${sessionid}" \
  "https://registry.jbei.org/rest/parts/${entryid}/traces"
