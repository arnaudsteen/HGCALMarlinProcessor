#!/bin/bash

source /Users/arnaudsteen/ilcsoft/init_ilcsoft.sh
export MARLIN_DLL=/Users/arnaudsteen/HGCAL/HGCALMarlinProcessor/lib/libhgcalMarlin.dylib

energy=$1
echo " .... Marlin process HGCAL display"

cat > LCIO.xml <<EOF
<marlin>
 <execute>
  <processor name="FSRStudyProcessor"/>
 </execute>

 <global>
  <parameter name="LCIOInputFiles">
   /data/lcio/calorimeterhit/multi_fsr_gamma${energy}GeV.slcio
  </parameter>
  <!-- limit the number of processed records (run+evt): -->  
  <!--parameter name="MaxRecordNumber" value="1"/-->  
  <!--parameter name="SkipNEvents" value="18000" /-->  
  <parameter name="SupressCheck" value="false" />  
  <!--parameter name="GearXMLFile"> gear_ldc.xml </parameter-->  
  <parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE  </parameter> 
  <!--parameter name="RandomSeed" value="1234567890" /-->
 </global>

 <processor name="FSRStudyProcessor" type="fsrStudyProcessor">
  <!--Name of the CalorimeterHit collection-->
  <parameter name="CollectionName" type="string" lcioInType="CalorimeterHit"> HGCALCalorimeterHit </parameter>
  <!--Name of the root output file-->
  <parameter name="OutputName" type="string" > fsr_gamma${energy}GeV.root </parameter>
  <parameter name="Geometry::NLayers" type="int"> 40 </parameter>
  <parameter name="Geometry::NPixelsPerLayer" type="int"> 128 </parameter>
  <parameter name="Geometry::Geometry::PixelSize" type="int"> 10.0 </parameter>
  <parameter name="InteractionFinder::MinSize" type="int"> 2 </parameter>
  <parameter name="InteractionFinder::MinEnergy" type="int"> 0.001 </parameter>
  <parameter name="Hough::NThetas" type="int"> 50 </parameter>
  <parameter name="Hough::MinimumNBins" type="int"> 10 </parameter>
  <parameter name="Hough::UseAnalogEnergy" type="bool"> true </parameter>
  <parameter name="Hough::MaxEnergy" type="float"> 0.0005 </parameter>
  <parameter name="Hough::MaximumNumberOfNeighboursForMip" type="int"> 2 </parameter>
 </processor>

</marlin>

EOF

Marlin LCIO.xml
#rm LCIO.xml
