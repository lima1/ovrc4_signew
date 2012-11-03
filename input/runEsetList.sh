#!/bin/sh
cp `locate createEsetList.R | grep library/curated` . 

R --vanilla "--args patientselection2.config.hold  eset.scaled.rda eset.scaled.log"  < createEsetList.R
R --vanilla "--args patientselection2binary.config.hold eset.binary.scaled.rda eset.binary.scaled.log"  < createEsetList.R
R --vanilla "--args patientselection2validation.config.hold eset.validation.scaled.rda eset.validation.scaled.log"  < createEsetList.R
R --vanilla "--args patientselection2debulking.config.hold eset.debulking.scaled.rda eset.debulking.scaled.log"  < createEsetList.R
R --vanilla "--args patientselection2early.config.hold eset.early.scaled.rda eset.early.scaled.log"  < createEsetList.R
#R --vanilla "--args patientselection2all.config.hold eset.all.scaled.rda eset.all.scaled.log"  < createEsetList.R
#R --vanilla "--args patientselection2allos.config.hold eset.allos.scaled.rda eset.allos.scaled.log"  < createEsetList.R
