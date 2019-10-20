name=$1
echo "============================"
echo "Running Example $name"
echo "============================"


thearray1=(mini rdwalk vprdwalk prspeed vrace ad2D vad1D american)
thearray2=(vmini single double vrdwalk prdwalk vprspeed race simple pollutant vad2D ad1D vamerican)
flag=1
for((i=1;i<=8;++i))
do
if [ $name == ${thearray1[i-1]} ];then
  flag=2
  sleep 1
fi
done


for((i=1;i<=12;++i))
do
if [ $name == ${thearray2[i-1]} ];then
  flag=3
  sleep 1
fi
done
one=1

echo "============================"
echo "Creating the Control Flow Graph..."
echo "============================"
sleep 1
if [ $flag -eq $one ];then

   ./CFG/cfg < Custom/Custom_inputs/$name.program Synthesis/${name}1.txt Synthesis/${name}2.txt
else
   ./CFG/cfg < Inputs/$name.program Synthesis/${name}1.txt Synthesis/${name}2.txt
fi
echo "============================"
echo "The Control Flow Graph is done."
echo "============================"
sleep 1


echo "============================"
echo "Example $name: Running the RSM-synthesis algorithm"
red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`
echo "${red}$(tput setab 7)This might take several minutes and MATLAB may appear to be frozen...${reset}"
echo "============================"
sleep 1
if [ $flag -eq $one ];then
   cp Custom/Custom_inputs/${name}config.txt Synthesis/${name}config.txt
else
   cp Inputs/${name}config.txt Synthesis/${name}config.txt
fi
cd Synthesis
echo "$name" | matlab -nodesktop -nosplash -r "run start.m"
cd ..

echo "============================"
echo "Example $name: RSM-synthesis is done."
echo "============================"
sleep 1

two=2

if [ $flag -eq $one ];then
  cp Synthesis/${name}_log.txt Custom/Custom_outputs/${name}_log.txt
  echo "============================"
  echo "Example $name: Final Output Written to Custom/Custom_outputs/${name}_log.txt"
  echo "============================"
  sleep 1
elif [ $flag -eq $two ];then
  cp Synthesis/${name}_log.txt Outputs/Table2/${name}_log.txt
  echo "============================"
  echo "Example $name: Final Output Written to Outputs/Table2/${name}_log.txt"
  echo "============================"
  sleep 1
else
  cp Synthesis/${name}_log.txt Outputs/Table3/${name}_log.txt
  echo "============================"
  echo "Example $name: Final Output Written to Outputs/Table3/${name}_log.txt"
  echo "============================"
  sleep 1
fi

echo "============================"
echo "Example $name: Removing Temporary Files"
echo "============================"
sleep 1
rm Synthesis/${name}1.txt > /dev/null 2>/dev/null
rm Synthesis/${name}2.txt > /dev/null 2>/dev/null
rm Synthesis/${name}config.txt > /dev/null 2>/dev/null
rm Synthesis/${name}_log.txt > /dev/null 2>/dev/null


echo "==========================="
echo "Example $name: Experiment Done"
echo "==========================="
sleep 1

sleep 1

