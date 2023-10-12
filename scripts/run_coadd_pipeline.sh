DATA_COLL=HSC/runs/RC2/w_2023_32/DM-40356
INJECT_COLL=u/pferguso/my-first-injection-catalog
REPO=/repo/main
OUTPUT=u/pferguso/my-first-stars
pipetask --long-log \
run --register-dataset-types \
--instrument lsst.obs.subaru.HyperSuprimeCam \
-b $REPO \
-i $DATA_COLL,$INJECT_COLL \
-o $OUTPUT \
-p DRP-RC2+injection_by_peter.yaml#step3 \
-d "instrument='HSC' AND skymap='hsc_rings_v1' AND
tract=9615 AND patch=3 AND band='i'"