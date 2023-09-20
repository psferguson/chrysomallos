$SOURCE_INJECTION_DIR/bin/generate_injection_catalog \
	-a 216.38 216.58  \
	-d -0.0960  0.0952 \
	-n 300 \
	-p source_type DeltaFunction \
	-p mag 13 15 17 19 21 \
    -f omg_coords_are_hard.csv \
    -b /repo/main \
    -i g r i \
    -o u/$USER/my-first-injection-catalog
