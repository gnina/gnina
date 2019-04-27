
#todo: add fixtures/depends clauses so these work properly when run in parallel
#https://crascit.com/2016/10/18/test-fixtures-with-cmake-ctest/

add_test(NAME grid0gpu COMMAND gninagrid -r files/rec.pdb -l files/lig.sdf -o time --gpu --time WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME grid1cpu COMMAND gninagrid -r files/rec.pdb -l files/lig.sdf -o lrcpu WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME grid1gpu COMMAND gninagrid -r files/rec.pdb -l files/lig.sdf -o lrgpu --gpu WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME grid1cmp COMMAND ./compare_bin.py lrcpu_0.48.28.binmap lrgpu_0.48.28.binmap WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

#random should get different result
add_test(NAME grid1rand COMMAND gninagrid -r files/rec.pdb -l files/lig.sdf -o lrrand --gpu  --random_rotation --random_translation 2 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME grid1randcmp COMMAND ./compare_bin.py lrgpu_0.48.28.binmap lrrand_0.48.28.binmap WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(grid1randcmp PROPERTIES WILL_FAIL TRUE)

add_test(NAME grid2gpu COMMAND gninagrid -r files/CC.xyz -l files/CC.xyz -o ccsmall --dimension 8 --resolution 0.25 --gpu --recmap files/recmap --ligmap files/ligmap WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME grid2cmp COMMAND ./compare_bin.py ccsmall_0.33.35.binmap files/ccsmall_0.33.35.binmap WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME gridbin COMMAND gninagrid -r files/CC.xyz -l files/CC.xyz -o ccbin --recmap files/recmap --ligmap files/ligmap WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME gridbincmp COMMAND ./compare_bin.py ccbin_0.48.35.binmap files/cc_0.48.35.binmap WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME griddx COMMAND gninagrid -r files/CC.xyz -l files/CC.xyz -o ccdx --dx --recmap files/recmap --ligmap files/ligmap WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME griddxcmpr COMMAND ./compare_dx.py ccdx_0_rec_AliphaticCarbonXSHydrophobe.dx files/ccdx_0_rec_AliphaticCarbonXSHydrophobe.dx WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME griddxcmpl COMMAND ./compare_dx.py ccdx_0_lig_AliphaticCarbonXSHydrophobe.dx files/ccdx_0_lig_AliphaticCarbonXSHydrophobe.dx WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME gridmap COMMAND gninagrid -r files/CC.xyz -l files/CC.xyz -o ccmap --gpu --map --recmap files/recmap --ligmap files/ligmap WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME gridmapcmpr COMMAND ./compare_map.py ccmap_0_rec_AliphaticCarbonXSHydrophobe.map files/ccmap_0_rec_AliphaticCarbonXSHydrophobe.map WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME gridmapcmpl COMMAND ./compare_map.py ccmap_0_lig_AliphaticCarbonXSHydrophobe.map files/ccmap_0_lig_AliphaticCarbonXSHydrophobe.map WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME gridbinary COMMAND gninagrid -r files/CC.xyz -l files/CC.xyz -o ccbin --binary_occupancy --gpu --dx --dimension 8 --recmap files/recmap --ligmap files/ligmap WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME gridbinarycmpr COMMAND ./compare_dx.py ccbin_0_rec_AliphaticCarbonXSHydrophobe.dx files/ccbin_0_rec_AliphaticCarbonXSHydrophobe.dx WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME gridbinarycmpl COMMAND ./compare_dx.py ccbin_0_lig_AliphaticCarbonXSHydrophobe.dx files/ccbin_0_lig_AliphaticCarbonXSHydrophobe.dx WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME griduser COMMAND gninagrid  -r files/C.xyz -l files/C.xyz -o ccgrid -g files/usergrid.dx  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME gridusercmp COMMAND ./compare_bin.py ccgrid_0.25.29.binmap files/ccgrid_0.25.29.binmap  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME gridsep COMMAND gninagrid -r files/A.xyz -l files/A.xyz -o ccsep --example_grid files/usergrid.dx  --separate WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME gridsepcmp COMMAND  ./compare_bin.py ccsep.25.14.binmap ccsep_0.25.14.binmap WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME gridsepnotcenter COMMAND bash -c "[ `od -f -w4 ccsep.25.14.binmap -v -Ad | grep 0031248 | awk '$2 < 0.5 {print \"done\"}'` == \"done\" ]"  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME gridcleanup COMMAND sh -c "rm *.binmap *.dx *.map" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
