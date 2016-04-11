#!/bin/bash

./HPGe_simulation gamma_hist_int_dl07.mac
sleep 1
./HPGe_simulation gamma_hist_int_dl08.mac
sleep 1
./HPGe_simulation gamma_hist_int_dl09.mac
sleep 1
./HPGe_simulation gamma_hist_int_dl10.mac
sleep 1
./HPGe_simulation gamma_hist_int_dl11.mac
sleep 1
./HPGe_simulation gamma_hist_int_dl12.mac
sleep 1
./HPGe_simulation gamma_hist_int_dl13.mac
sleep 1
./HPGe_simulation gamma_hist_int_dl14.mac
sleep 1
./HPGe_simulation gamma_hist_int_dl15.mac
sleep 1
./HPGe_simulation gamma_hist_int_dl16.mac
sleep 1
./HPGe_simulation gamma_hist_int_dl17.mac
sleep 1
./HPGe_simulation gamma_hist_int_dl18.mac
sleep 1
mv  ./pointd90mmdl*.root ~/Nutstore/ROOTScript/DataAnalysis/data/pointd90dl_new

exit 1
