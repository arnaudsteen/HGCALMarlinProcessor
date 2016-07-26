#!/bin/bash

source /Users/arnaudsteen/ilcsoft/init_ilcsoft.sh
export MARLIN_DLL=/Users/arnaudsteen/HGCAL/HGCALMarlinProcessor/lib/libhgcalMarlin.dylib


nparticle=$1
seed=$2
echo " .... Marlin process HGCAL display"

cat > LCIO.xml <<EOF
<marlin>
 <execute>
  <processor name="HGCALMuonFinder"/>
 </execute>

 <global>
  <parameter name="LCIOInputFiles">
   /data/lcio/calorimeterhit/muonPlusJet${nparticle}_${seed}.slcio
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
  <parameter name="OutputName" type="string" > muonPlusJet${nparticle}_${seed}.root </parameter>
  <!--Number of steps to loop over theta values-->
  <parameter name="LayerZPosition" type="IntVec"> -209.55 -200.85 -187.75 -179.05 -165.95 -157.25 -144.15 -135.45 -122.35 -113.65 -100.55 -91.85 -78.75 -70.05 -56.95 -48.25 -35.15 -26.45 -13.35 -4.65 8.45 17.15 30.25 38.95 52.05 60.75 73.85 82.55 95.65 104.35 117.45 126.15 139.25 147.95 161.05 169.75 182.85 191.55 204.65 213.35 </parameter>
  <parameter name="Geometry::NLayers" type="int"> 40 </parameter>
  <parameter name="Geometry::NPixelsPerLayer" type="int"> 192 </parameter>
  <parameter name="Geometry::Geometry::PixelSize" type="int"> 10.0 </parameter>
  <parameter name="Geometry::Geometry::FirstLayerZ" type="float"> -209.55 </parameter>
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
