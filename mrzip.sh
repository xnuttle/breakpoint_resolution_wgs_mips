for i in $(ls|grep 'qseq');
  do echo "gzip $i &"
done >> zip_commands.sh
echo "wait" >> zip_commands.sh
echo "rm -f zip_commands.sh" >> zip_commands.sh
