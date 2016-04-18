#!/bin/bash

source /Users/arnaudsteen/ilcsoft/init_ilcsoft.sh
export MARLIN_DLL=/Users/arnaudsteen/HGCAL/HGCALMarlinProcessor/lib/libhgcalMarlin.dylib

particle=$1
energy=$2
seed=$3
echo " .... Marlin process HGCAL display"

cat > LCIO.xml <<EOF
<marlin>
 <execute>
  <processor name="ShowerProcessor"/>
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

 <processor name="ShowerProcessor" type="hgcalShowerProcessor">
  <!--Name of the CalorimeterHit collection-->
  <parameter name="CollectionName" type="string" lcioInType="CalorimeterHit"> HGCALCalorimeterHit </parameter>
  <!--Name of the root output file-->
  <parameter name="OutputName" type="string" > single_${particle}_${energy}GeV_${seed}.root </parameter>
 </processor>

</marlin>

EOF

Marlin LCIO.xml
#rm LCIO.xml
