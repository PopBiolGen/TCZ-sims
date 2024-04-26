setwd("~/final_natural_water_files/PBS_Files")

script.file<-'~/evo-dispersal/art_wbdies/VRD/spread_model.R'

for (ss in 1:100) {
  
    
    	fid<-ss
	    ##create the sh file
	    zz = file(paste("VRD_ABC", fid,'.sh',sep=''),'w')
  	  cat('##################################\n',file=zz)
	    cat('#!/bin/sh\n',file=zz)
	    cat('cd $PBS_O_WORKDIR\n',file=zz)
  	  cat("R CMD BATCH --no-save --no-restore '--args File.ID=", fid, "' ", sep="", file=zz)
	    cat(script.file, " ", paste("VRD_ABC", fid,'.Rout',sep=''), "\n", sep="",file=zz)
	    cat('##################################\n',file=zz)
  	  close(zz)
			
    	#submit the job
	    system(paste("qsub -m n VRD_ABC", fid,".sh",sep=""))
  
}
