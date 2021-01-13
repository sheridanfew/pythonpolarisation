CWD=$(pwd)
for file in */*coul*/*30*py
do
  cd $(dirname $file)
  python *py
  cd $CWD
done

for file in */*coul*/*30*csv
do
  cat $file >> coulombints.csv
  echo -en "\n" >> coulombints.csv
done

