ROSETTA_BIN="/home/liusai/software/rosetta/rosetta.binary.ubuntu.release-371/main/source/bin"
INPUT_DIR="./data"
WORK_DIR="./results_full"
WEIGHTS="ref2015"
INTERFACE="A_B"
NSTRUCT=1
NPROC=31 
mkdir -p "$WORK_DIR"
process_one() {
    pdb="$1"
    filename=$(basename "$pdb" .pdb)
    outdir="$WORK_DIR/$filename"
    mkdir -p "$outdir"

    echo "🔧 Relaxing: $filename"
    "$ROSETTA_BIN/relax.static.linuxgccrelease" \
        -s "$pdb" \
        -nstruct $NSTRUCT \
        -relax:fast \
        -relax:constrain_relax_to_start_coords \
        -relax:coord_cst_stdev 0.3 \
        -relax:coord_constrain_sidechains \
        -score:weights "$WEIGHTS" \
        -out:pdb \
        -out:path:all "$outdir" \
        -out:file:scorefile "$outdir/relax_score.sc" \
        -ignore_unrecognized_res \
        -overwrite

    relaxed_pdb="$outdir/${filename}_0001.pdb"
    echo " Analyzing interface: $filename"
    "$ROSETTA_BIN/InterfaceAnalyzer.static.linuxgccrelease" \
        -s "$relaxed_pdb" \
        -interface $INTERFACE \
        -score:weights "$WEIGHTS" \
        -pack_separated \
        -compute_packstat \
        -ignore_unrecognized_res \
        -out:file:scorefile "$outdir/interface_score.sc" \
        -overwrite
}
export -f process_one
export ROSETTA_BIN WEIGHTS WORK_DIR NSTRUCT INTERFACE
