label=$1
array1=(mini rdwalk vprdwalk prspeed vrace ad2D vad1D american)

array2=(vmini single double vrdwalk prdwalk vprspeed race simple pollutant vad2D ad1D vamerican)

num1=${#array1[*]}
num2=${#array2[*]}

if [ $label == "table2" ];then
  echo "================================="
  echo "Data generated for Table 2 in Section 9 (see the paper)"
  echo "================================="
  for((i=1;i<=$num1;++i))
  do	
     ./run-example.sh ${array1[i-1]}
  done
elif [ $label == "table3" ];then
  echo "================================="
  echo "Data generated for Table 3 in Section 9 (see the paper)"
  echo "================================="
  for((i=1;i<=$num2;++i))
  do	
     ./run-example.sh ${array2[i-1]}
  done
else
  echo "================================="
  echo "Data generated for Table 2 and Table 3 in Section 9 (see the paper)"
  echo "================================="
  for((i=1;i<=$num1;++i))
  do	
     ./run-example.sh ${array1[i-1]}
  done

  for((i=1;i<=$num2;++i))
  do	
     ./run-example.sh ${array2[i-1]}
  done
fi

