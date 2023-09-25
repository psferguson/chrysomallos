DATA_COLL=HSC/runs/RC2/w_2023_32/DM-40356
INJECT_COLL=u/pferguso/full_dwarf_test/1_lim_32
REPO=/repo/main
OUTPUT=u/pferguso/full_dwarf_test/1_lim_32/step3
pipetask --long-log \
run --register-dataset-types \
--instrument lsst.obs.subaru.HyperSuprimeCam \
-b $REPO \
-i $DATA_COLL,$INJECT_COLL \
-o $OUTPUT \
-p pipelines/DRP-RC2+injection_by_peter.yaml#step3 \
-d "instrument='HSC' AND skymap='hsc_rings_v1' AND
tract=9615 AND patch=3" \
-j 3


DATA_COLL=HSC/runs/RC2/w_2023_32/DM-40356
INJECT_COLL=u/pferguso/full_dwarf_test/1_lim_28
REPO=/repo/main
OUTPUT=u/pferguso/full_dwarf_test/1_lim_28/step3
pipetask --long-log \
run --register-dataset-types \
--instrument lsst.obs.subaru.HyperSuprimeCam \
-b $REPO \
-i $DATA_COLL,$INJECT_COLL \
-o $OUTPUT \
-p pipelines/DRP-RC2+injection_by_peter.yaml#step3 \
-d "instrument='HSC' AND skymap='hsc_rings_v1' AND
tract=9615 AND patch=3" \
-j 3