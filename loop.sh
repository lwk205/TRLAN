for x in `ls | grep th_`; do
echo $x
cat $x | grep TIME
done

