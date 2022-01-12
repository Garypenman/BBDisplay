#!/bin/csh
echo " setting DATA_DIR etc"

export SBS=$HOME/sbs/sbs_devel/install
export SBSOFFLINE=$HOME/sbs/sbs_devel/install
export SBS_REPLAY=$HOME/sbs/sbs_devel/SBS-replay

export DISPLAY_DIR=$HOME/sbs/BBDisplay
export DATA_DIR=$HOME/sbs/data

export DB_DIR=$SBS_REPLAY/DB
#export ANALYZER_CONFIGPATH=$SBS_REPLAY/replay
export ANALYZER_CONFIGPATH=$DISPLAY_DIR/replay

export OUT_DIR=$DISPLAY_DIR/Rootfiles
export LOG_DIR=$DISPLAY_DIR/AnalysisLogs

export LD_LIBRARY_PATH=$SBS/lib64:$LD_LIBRARY_PATH
export PATH=$SBS/bin:$PATH
