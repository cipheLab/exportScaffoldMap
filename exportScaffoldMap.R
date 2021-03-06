library("igraph")
library("scaffold")
library("stringr")

library("utils")
library("flowCore")
library("Biobase")

### demo

cytofCore.updateFlowFrameKeywords = function(flowFrame){
  
  row.names(flowFrame@parameters) <- paste0("$P",c(1:length(row.names(flowFrame@parameters))))
  params = parameters(flowFrame)
  pdata = pData(params)
  for (i in 1:ncol(flowFrame)){
  	n = paste0("$P",i,"N");
    s = paste0("$P",i,"S");
    r = paste0("$P",i,"R");
    b = paste0("$P",i,"B");
    e = paste0("$P",i,"E");
    keyval=list();
    label <- pData(flowFrame@parameters)[,"desc"][i]
    if(is.na(label)) {label <- colnames(flowFrame)[i] }
    keyval[[n]] = colnames(flowFrame)[i] 
    keyval[[s]] = label       
    keyval[[r]] = ceiling(max(exprs(flowFrame)[,i]))
    keyval[[b]] = 32
    keyval[[e]] = "0,0"
    keyword(flowFrame) = keyval;
    
    pdata[i,"minRange"]=min(exprs(flowFrame)[,i])
    pdata[i,"maxRange"]=max(exprs(flowFrame)[,i])
    
  }
  pData(params)=pdata
  parameters(flowFrame)=params
  row.names(flowFrame@parameters) <- paste0("$P",c(1:length(row.names(flowFrame@parameters))))
  # keyval[["$DATATYPE"]] <- "F"
  return(flowFrame)
}

exportScaffoldMap <- function(inputPath=NULL, outputPath=NULL, trans.value = 5,
															new.fcs=TRUE, 
															new.txt=FALSE, 
															pop.csv=TRUE,
															cluster.intensity=FALSE){

	if(is.null(outputPath)){outputPath<-inputPath}

	scaffold.Data <- list.files(inputPath, pattern=".scaffold", full.names=TRUE, recursive=FALSE)
	list.txt <- list.files(inputPath, pattern=".txt$", full.names=TRUE, recursive=FALSE)
	list.fcs <- list.files(inputPath, pattern=".fcs$", recursive=FALSE, full.names=TRUE)
	list.data <- list.files(inputPath, pattern=".RData$", recursive=FALSE, full.names=TRUE)

	if(length(scaffold.Data)>1){
		print("Warning : more than one scaffold file found, just first used")
		scaffold.Data <- scaffold.Data[1] 
	}

	data <- scaffold:::my_load(scaffold.Data)
	all.table <- lapply(c(1:length(data$graphs)), function(i){
		G <- data$graphs[[i]]
		x <- V(G)$x
		y <- V(G)$y
		popsize <- V(G)$popsize
		if(is.null(popsize)){popsize <- (rep(1,length(x)))}
		label <- V(G)$groups
		
		tab <- data.frame(x,y,popsize,label)
		return(tab)
	})

	node.index.x <- which(V(data$graphs[[1]])$type==1)

	table.node <- data.frame(all.table[[1]][node.index.x,c(1,2)])
	G <- data$graphs[[1]]
	pop.name <- unlist(lapply(c(1:dim(table.node)[1]),function(x){return(names(G[[x]]))}))
	table.node <- cbind(table.node, pop.name, c(1:nrow(table.node)))
	colnames(table.node) <- c("Scaffold.x","Scaffold.y","pop","popID.Scaffold")
	write.csv(table.node, paste0(outputPath,"node_landmark.coordinate.csv"), row.names=FALSE)

	print("Coordinate csv in progress...")
	lapply(c(1:length(list.txt)), function(x){

		name <- basename(list.txt[[x]])
		id <- which(names(data$graphs)==name)
		if(length(id)==0){
			return(NULL)
		} else {
			d <- all.table[[id]]

			txt <- read.csv(list.txt[[x]],header=TRUE,sep="\t",check.names = FALSE)
			txt <- txt[,-(which(colnames(txt)%in%c("cellType","sample","popsize")))]
			tot <- dim(read.FCS(list.fcs[x]))[1]

			
			table <- data.frame(d[-node.index.x,])
			temp <- as_data_frame(data$graphs[[id]])
			pop <- temp[which(temp[,"edge_type"]=="highest_scoring"),"from"]

			table <- cbind(table,table[,3],pop,txt)
			table[,3] <- (table[,3]/tot)*100
			colnames(table) <- c("Scaffold.x","Scaffold.y","PercentOfTotal","clusterID.Scaffold","Cluster.Size","pop",colnames(txt)) #clusterID,x,y,Events,Percentile,MFIs,
			write.csv(table, paste0(outputPath,sub(".fcs","",basename(list.fcs[x])),".coordinate.csv"), row.names=FALSE)
			print(paste0("Generate CSV : ",x))
		}
	})
	print("Coordinate csv Done")

	if(new.fcs == TRUE){
		print("Generate new FCS...")
		landmark <- read.csv(paste0(outputPath,"node_landmark.coordinate.csv"))
		lapply(c(1:length(list.fcs)),function(i){
			fcs <- read.FCS(list.fcs[i])
			name <- basename(list.txt[[i]])
			id <- which(names(data$graphs)==name)
			if(length(id)==0){
				return(NULL)
			} else {
				d <- all.table[[id]]
				table <- read.csv(paste0(outputPath,sub(".fcs","",basename(list.fcs[i])),".coordinate.csv"))
				rdata <- scaffold:::my_load(list.data[i])
				new_col.1 <- matrix(rdata[,"cellType"], nrow = nrow(rdata), ncol = 1, dimnames = list(NULL, "clusterID"))
				
				new_col.2 <- as.vector(table[new_col.1,"pop"])
				pop.index <- unlist(lapply(new_col.2,function(j){return(landmark[which(j==landmark[,"pop"]),"popID.Scaffold"])}))
				tmp <- matrix(landmark[pop.index,"popID.Scaffold"], ncol = 1, dimnames = list(NULL, "popID.Scaffold"))

				fcs <- flowCore::cbind2(fcs, new_col.1)
				fcs <- flowCore::cbind2(fcs, tmp)
				row.names(pData(fcs@parameters)) <- paste0("$P",c(1:length(row.names(fcs@parameters))))
				fcs@description[[paste0("$P",(dim(fcs)[2]-1),"R")]] <- 262144
				fcs@description[[paste0("$P",(dim(fcs)[2]),"R")]] <- 262144
				# fcs <- cytofCore.updateFlowFrameKeywords(fcs)	
				print(paste0("Generate new FCS : ",i))
				write.FCS(fcs, paste0(outputPath,sub(".fcs","",basename(list.fcs[i])),"_celltype.fcs"),delimiter="#")
			}
		})
		print("Generate new FCS Done")
	}

	if(new.txt == TRUE){
		print("Generate txt file...")
		lapply(c(1:length(list.fcs)),function(x){

			fcs <- read.FCS(list.fcs[x])
			name <- basename(list.txt[[x]])
			id <- which(names(data$graphs)==name)
			if(length(id)==0){
				return(NULL)
			} else {

				print(paste("Generate txt file : ",x))
				fcs <- read.FCS(paste0(outputPath,sub(".fcs","",basename(list.fcs[x])),"_celltype.fcs"),emptyValue=TRUE)
				asinhTrans <- arcsinhTransform(transformationId="ln-transformation", a=trans.value, b=1, c=1)
				param <- colnames(fcs)[0:(dim(fcs)[2]-4)]
				translit <- transformList(param, asinhTrans)
				fcs <- transform(fcs, translit)
				mat <- fcs@exprs
				write.csv(mat, paste0(outputPath,sub(".fcs","",basename(list.fcs[x])),"_transform_cellType.csv"), row.names=FALSE)
			}
		})
		print("Generate txt DONE.")
	}

	if(pop.csv == TRUE){
		print("Generate csv pop file...")
		lapply(c(1:length(list.fcs)),function(x){

			fcs <- read.FCS(list.fcs[x])
			name <- basename(list.txt[[x]])
			id <- which(names(data$graphs)==name)
			if(length(id)==0){
				return(NULL)
			} else {

				print(paste0("Generate csv pop file : ",x))
				fcs <- read.FCS(paste0(outputPath,sub(".fcs","",basename(list.fcs[x])),"_celltype.fcs"),emptyValue=FALSE)
				rdata <- scaffold:::my_load(list.data[x])
				res <- lapply(unique(fcs@exprs[,"popID.Scaffold"]),function(y){
				  print(y)
					mat <- fcs@exprs[which(fcs@exprs[,"popID.Scaffold"]==y),]
					if(dim(mat)[1]< 2 || is.null(dim(mat))){return(c(0,0))}
					return(c(apply(mat,2,mean),dim(mat)[1]))
				})
				table <- do.call(rbind, res)
				colnames(table)[dim(table)[2]] <- "Count"
				write.csv(table, paste0(outputPath,sub(".fcs","",basename(list.fcs[x])),"_popMFI.csv"))
			}
		})
		print("Generate csv pop DONE.")
	}

	if(cluster.intensity == TRUE){
		print("Generate intensity csv file...")
		lapply(c(1:length(list.fcs)),function(x){
			print(paste0("Generate intensity csv file : ",x))
			fcs <- read.FCS(paste0(outputPath,sub(".fcs","",basename(list.fcs[x])),"_celltype.fcs"),emptyValue=FALSE)
			res <- lapply(unique(fcs@exprs[,"clusterID"]),function(y){
				mat <- fcs@exprs[which(fcs@exprs[,"clusterID"]==y),]
				if(dim(mat)[1]< 2 || is.null(dim(mat))){return(c(0,0))}
				return(c(apply(mat,2,median),dim(mat)[1]))
			})
			table <- do.call(rbind, res)
			colnames(table)[dim(table)[2]] <- "Count"
			write.csv(table, paste0(outputPath,sub(".fcs","",basename(list.fcs[x])),"_clusterIntensity.csv"))
		})
		print("Generate intensity cluster DONE.")
	}
}


exportScaffoldMap.run <- function(inputPath=dirname(file.choose()),new.fcs=TRUE,new.txt=FALSE, pop.csv=TRUE, cluster.intensity = FALSE){
	exportScaffoldMap(paste0(normalizePath(inputPath,winslash = "/"),"/"),new.fcs=new.fcs,new.txt=new.txt,pop.csv=pop.csv,cluster.intensity=cluster.intensity)
}
