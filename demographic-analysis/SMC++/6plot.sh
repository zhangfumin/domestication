#date_200831_plot
#wangmeixia

export PATH=/home/gesong/wangmeixia/datapool/software/miniconda3/bin:$PATH 
source activate /home/gesong/wangmeixia/datapool/software/miniconda3/
conda activate smcpp

smc++ plot all_estimate.pdf    ./ind1/model_time/model.final.json ./aus/model_time/model.final.json ./aro/model_time/model.final.json ./ray/model_time/model.final.json ./trj3/model_time/model.final.json ./tej/model_knots/model.final.json ./NIV1/model_time/model.final.json ./NIV2/model_time/model.final.json ./RUF1/model_time/model.final.json ./RUF2_1/model_time/model.final.json -c -x 50  20000 -y 300  300000
smc++ plot sativa_estimate.pdf ./ind1/model_time/model.final.json ./aus/model_time/model.final.json ./aro/model_time/model.final.json ./ray/model_time/model.final.json ./trj3/model_time/model.final.json ./tej/model_knots/model.final.json -c -x 50  20000 -y 300  300000
smc++ plot wild_estimate.pdf    ./NIV1/model_time/model.final.json ./NIV2/model_time/model.final.json ./RUF1/model_time/model.final.json ./RUF2_1/model_time/model.final.json -c -x 50  20000 -y 300  300000
