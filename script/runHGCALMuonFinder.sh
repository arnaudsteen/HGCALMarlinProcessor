#!/bin/bash

source /Users/arnaudsteen/ilcsoft/init_ilcsoft.sh
export MARLIN_DLL=/Users/arnaudsteen/HGCAL/HGCALMarlinProcessor/lib/libhgcalMarlin.dylib

seed=$1
echo " .... Marlin process HGCAL display"

cat > LCIO.xml <<EOF
<marlin>
 <execute>
  <processor name="HGCALMuonFinder"/>
 </execute>

 <global>
  <parameter name="LCIOInputFiles">
   /data/lcio/calorimeterhit/mu-_jet_${seed}.slcio
  </parameter>
  <!-- limit the number of processed records (run+evt): -->  
  <!--parameter name="MaxRecordNumber" value="1000"/-->  
  <!--parameter name="SkipNEvents" value="18000" /-->  
  <parameter name="SupressCheck" value="false" />  
  <!--parameter name="GearXMLFile"> gear_ldc.xml </parameter-->  
  <parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE  </parameter> 
  <!--parameter name="RandomSeed" value="1234567890" /-->
 </global>

 <processor name="HGCALMuonFinder" type="hgcalMuonFinder">
  <!--hgcalDisplayProcessor calculates a HGCAL display and pad multiplcity for each SHDCAL layer-->
  <!--Name of the CalorimeterHit collection-->
  <parameter name="CollectionName" type="string" lcioInType="CalorimeterHit"> HGCALCalorimeterHit </parameter>
  <!--Name of the root output file-->
  <parameter name="OutputName" type="string" > mu-_jet_${seed}.root </parameter>
  <!--Number of steps to loop over theta values-->
  <parameter name="Geometry::NLayers" type="int"> 40 </parameter>
  <parameter name="Geometry::NPixelsPerLayer" type="int"> 192 </parameter>
  <parameter name="Geometry::Geometry::PixelSize" type="int"> 10.0 </parameter>
  <parameter name="Geometry::Geometry::FirstLayerZ" type="float"> -209.60 </parameter>
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
rm LCIO.xml
