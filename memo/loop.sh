for x in `ls | grep res`; do
echo $x
cat $x | grep TIME
done

