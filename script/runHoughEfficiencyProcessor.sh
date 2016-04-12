#!/bin/bash

source /Users/arnaudsteen/ilcsoft/init_ilcsoft.sh
export MARLIN_DLL=/Users/arnaudsteen/HGCAL/HGCALMarlinProcessors/lib/libhgcalMarlin.dylib

particle=$1
energy=$2
seed=$3
echo " .... Marlin process HGCAL display"

cat > LCIO.xml <<EOF
<marlin>
 <execute>
  <processor name="HoughEfficiency"/>
 </execute>

 <global>
  <parameter name="LCIOInputFiles">
   /data/lcio/calorimeterhit/single_${particle}_${energy}GeV_${seed}.slcio
  </parameter>
  <!-- limit the number of processed records (run+evt): -->  
  <!--parameter name="MaxRecordNumber" value="1000"/-->  
  <!--parameter name="SkipNEvents" value="18000" /-->  
  <parameter name="SupressCheck" value="false" />  
  <!--parameter name="GearXMLFile"> gear_ldc.xml </parameter-->  
  <parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE  </parameter> 
  <!--parameter name="RandomSeed" value="1234567890" /-->
 </global>

 <processor name="HoughEfficiency" type="houghEfficiencyProcessor">
  <!--hgcalDisplayProcessor calculates a HGCAL display and pad multiplcity for each SHDCAL layer-->
  <!--Name of the CalorimeterHit collection-->
  <parameter name="CollectionName" type="string" lcioInType="CalorimeterHit"> HGCALCalorimeterHit </parameter>
  <!--Name of the root output file-->
  <parameter name="OutputName" type="string" > single_${particle}_${energy}GeV_${seed}.root </parameter>
  <!--Number of steps to loop over theta values-->
  <parameter name="Hough::NThetas" type="int"> 50 </parameter>
 </processor>

</marlin>

EOF

Marlin LCIO.xml
#rm LCIO.xml
