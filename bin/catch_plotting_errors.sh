#!/bin/bash -e

if [ -s $1 ]; then
    echo "Caught a plotting error"

    NAME="Consensus quality"
    PRIORITY=91
    ERROR=`(cat $1)`
    PAGE_CONTENTS="<h1>Plotting for $NAME failed</h1><b>ERROR MESSAGE</b><br><p style="font-family:monospace"><FONT color="#DC381F">$ERROR</FONT></p>"
    STATUS=1

    JSON_FMT='{"%s": {"name":"%s","priority":%s,"script":"%s","div":"%s"},"additional_info": {"exit_status":%s}}\n'
    printf "$JSON_FMT" "$NAME" "$NAME" "$PRIORITY" "$PAGE_CONTENTS" "$PAGE_CONTENTS" "$STATUS" > $2

else
    echo "No plotting errors"
fi