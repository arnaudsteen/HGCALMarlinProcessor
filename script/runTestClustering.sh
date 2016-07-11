#!/bin/bash

source /Users/arnaudsteen/ilcsoft/init_ilcsoft.sh
export MARLIN_DLL=/Users/arnaudsteen/HGCAL/HGCALMarlinProcessor/lib/libhgcalMarlin.dylib

echo " .... Marlin process HGCAL display"

path="/data/lcio/calorimeterhit"

particle="mu-"
energy=10
seed=0
prefix="muon"

if test $# -eq 3
then 
    particle=$1
    energy=$2
    seed=$3
    file="single_${particle}_${energy}GeV_${seed}.slcio"
    prefix="${particle}"
else test $# -eq 1
     energy=$1
     file="multi_fsr_gamma${energy}GeV.slcio"
     prefix="multi_fsr_gamma${energy}GeV"
fi

list=`ls $path/$file`
echo "list : "$list

cat > LCIO.xml <<EOF
<marlin>
 <execute>
  <processor name="testClusteringProcessor"/>
 </execute>
 
 <global>
  <parameter name="LCIOInputFiles">
   ${path}/${file}
  </parameter>
  <!-- limit the number of processed records (run+evt): -->  
  <parameter name="MaxRecordNumber" value="50"/>  
  <!--parameter name="SkipNEvents" value="18000" /-->  
  <parameter name="SupressCheck" value="false" />  
  <!--parameter name="GearXMLFile"> gear_ldc.xml </parameter-->  
  <parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE  </parameter> 
  <!--parameter name="RandomSeed" value="1234567890" /-->
 </global>
 
 <processor name="testClusteringProcessor" type="testClusteringProcessor">
  <!--testClusteringProcessor displays events in HGCAL-->
  <!--Name of the CalorimeterHit collection-->
  <parameter name="CollectionName" type="string" lcioInType="CalorimeterHit"> HGCALCalorimeterHit </parameter>
  <!--Name of the root output file-->
  <parameter name="PrefixPlotName" type="int"> ${prefix} </parameter>
  <parameter name="PauseAfterDraw" type="bool"> true </parameter>
  <parameter name="Geometry::Geometry::FirstLayerZ" type="float"> -144.15 </parameter>
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
  <parameter name="Cluster3D::MaxLongitudinal" type="int"> 2 </parameter>
  <parameter name="Cluster3D::MaxTransverseDistance" type="float"> 100.0 </parameter>
 </processor>
</marlin>

EOF

Marlin LCIO.xml
#rm LCIO.xml
