for F in $( ls ./*.sh ); do echo $F && shellcheck.sh $F; done
