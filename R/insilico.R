#' @importFrom data.table fread
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom plyr ldply


degenerate.pairs <- list(
	"A" = "A",
	"C" = "C",
	"G" = "G",
	"T" = "T",
	"R" = c("A", "G"),
	"Y" = c("C", "T"),
	"M" = c("A", "C"),
	"K" = c("G", "T"),
	"S" = c("C", "G"),
	"W" = c("A", "T"),
	"H" = c("A", "C", "T"),
	"B" = c("C", "G", "T"),
	"V" = c("A", "C", "G"),
	"D" = c("A", "G", "T"),
	"N" = c("A", "C", "G", "T")
	);

demux <- function(seq) {
	seq <- unlist(strsplit(seq,""))
	g <- expand.grid(degenerate.pairs[match(seq, names(degenerate.pairs))])
	g <- unlist(apply(g,1,function(x){paste(x,collapse="")}))
	return(g)
	
	}

#' Create a BLAST Database
#'
#' Wrapper function for `makeblastdb`
#' @param x an object of class XStringSet
#' @param dbtype character vector indicating molecule type  (\code{"nucl"} or \code{"prot"})
#' @return String giving the path to the BLAST database
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' @export

makeblastdb <- function(x,dbtype="nucl") {
	fn <- tempfile();
	writeXStringSet(x,fn);
	ex <- paste(Sys.which("makeblastdb"),"-in",fn,"-dbtype",dbtype)
	system(ex)
	fn
	}

#' Create a BLAST Database for Insilico PCR primers
#' 
#' Demultiplexes degenerate bases in a set of forward and reverse primers and creates a blast database
#' @param x a vector or list of length 2 containing a forward primer (5'->3' on plus strand) and reverse primer (5'->3' on minus strand) 
#'
#' @return String giving the path to the BLAST database
#' @examples
#' primers <- c("GGCTGGATCACCTCCTT","GCCWAGGCATCCDCC")
#' db <- makeprimerdb(primers)
#'
#'
#'
#' @export

makeprimerdb <- function(x){
	names(x) <- c("for","rev");
	s <- DNAStringSet(unlist(lapply(x,demux)));
	makeblastdb(s)
	}
#' NCBI BLAST (Basic Local Alignment Search Tool)	
#' @param query either a string pointing to the query file, an object of class DNAbin, or an object of class XStringSet
#' @param db string pointing to the BLAST database, or if the switch `remote` is provided in `args`, the name of the remote NCBI database
#' @param args A named list containing the BLAST program arguments to be used in a key/value manner. There is no need to use leading dashes on the argument names. SEE EXAMPLES
#' @param type string with the name of the BLAST program to use. e.g. `blastp` or `blastn`
#' @param gz boolean. For gzip'd query files. if `TRUE`, the query file is uncompressed and read stdin into the BLAST program
#' @return An object of class `data.table`. Only format `6` (tabular format) for `outfmt` is supported.
#' @examples
#' ##Run local primer blast
#'
#' query <- readDNAStringSet("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/010/525/GCA_000010525.1_ASM1052v1/GCA_000010525.1_ASM1052v1_genomic.fna.gz")
#'
#' primers <- c("GGCTGGATCACCTCCTT","GCCWAGGCATCCDCC")
#' db <- makeprimerdb(primers)
#'
#' #Define BLAST arguments: these are the same that can be found in `RinsilicoPCR::blast.params$primers`
#' args <- list(
#' 		outfmt="'6 qacc sacc pident length qlen slen qstart sstart qend send evalue'",
#' 		max_target_seqs=20,
#' 		word_size=8
#' 		)
#'
#' blastres <- blast(query,db,args)
#' 
#'
#'
#' ##Run remote blast on NCBI blastdb `nt`
#'
#' query <- readDNAStringSet("query.fa")
#' 
#' #Define BLAST arguments: these are the same that can be found in `RinsilicoPCR::blast.params$remote`
#' args <- list(
#'		outfmt="'6 qacc sacc pident length qlen slen qstart sstart qend send evalue'",
#'		max_target_seqs=20,
#'		word_size=11,
#'		remote=""
#'		)
#'
#' blastres <- blast(query,db="nt",args)
#' 



#' @export

blast <- function(query,db,args,type="blastn",gz=FALSE){
	fn <- list(
		"character"=function(){
			query
			},
		"DNAStringSet"=function(){
			fn <- tempfile()
			writeXStringSet(query,fn)
			fn
			},
		"DNAbin"=function(){
			fn <- tempfile()
			write.dna(query,fn)
			fn
			}
		)[[class(query)[1]]]()
	
	read_method <- if(gz) "gzip -dc" else "cat"
	
	args_s <- paste(sapply(names(args),function(y){paste("-",y," ",args[[y]],sep="")}),collapse=" ")
	
	
	exec <- paste(read_method,fn,"|",Sys.which(type),"-db",db,args_s)
	structure(
		fread(exec),
		names=strsplit(gsub("[\\'\\\"]","",args$outfmt),"[ ][ ]*")[[1]][-1],
		query=fn
		)
	}
	
#' @export
blast.params <- list(
	primers = list(
		outfmt="'6 qacc sacc pident length qlen slen qstart sstart qend send evalue'",
		max_target_seqs=20,
		word_size=8
		),
	remote= list(
		outfmt="'6 qacc sacc pident length qlen slen qstart sstart qend send evalue'",
		max_target_seqs=20,
		word_size=11,
		remote=""
		)
	)
	

#' Find primer pairs in BLAST results
#' @param x data.frame of BLAST results
#' @param filt Logical expression used to filter mismatches from primer BLAST results. Defaults to `pident >= 80 & length/slen >= .8`, (percent identity > 80\% and match length > 80\% of primer).
#' @param size.range numeric vector of length 2 giving the range of expected amplicon length.
#' @return list of results for each query sequence containing:
#' @return dir : orientation of amplicon
#' @return blastres : BLAST table results for primer pair
#' @return location : start and end coordinates for the amplicon
#' @return size : size of the amplicon
#' @examples
#' ##Run local primer blast
#'
#' query <- readDNAStringSet("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/010/525/GCA_000010525.1_ASM1052v1/GCA_000010525.1_ASM1052v1_genomic.fna.gz")
#'
#' primers <- c("GGCTGGATCACCTCCTT","GCCWAGGCATCCDCC")
#' db <- makeprimerdb(primers)
#'
#' #Define BLAST arguments: these are the same that can be found in `RinsilicoPCR::blast.params$primers`
#' args <- list(
#' 		outfmt="'6 qacc sacc pident length qlen slen qstart sstart qend send evalue'",
#' 		max_target_seqs=20,
#' 		word_size=8
#' 		)
#'
#' blastres <- blast(query,db,args)
#' 
#' pairs <- primer.insilico(blastres)
#' @export

primer.insilico <- function(x,filt=(pident >= 80 & length/slen >= .8),size.range=c(50,1200)) {
	
	filt <- deparse(substitute(filt))
	f <- with(x,eval(parse(text=filt)))
	
	r <- x[f,]
	
	
	r$ptype <- c(1,-1)[match(gsub("[0-9]*$","",r$sacc),c("for","rev"))]
	r$dir <- sign(r$send - r$sstart)
	
	#apply for each sequence in fasta
	o <- lapply(split(r,r$qacc),function(x){
		
		#remove redundant blast hits to other degenerate primers
		x <- x[order(x$qstart),]
		d <- which(c(NA,diff(x$qstart)) <= 7)
		
		if(length(d)) {
			dup <- unlist(lapply(split(d,cumsum(c(2,diff(d)) != 1)),function(y) {
				y <- c(min(y)-1,y)
				y[order(x$evalue[y])][-1]
				}),use.names=FALSE)
			x <- x[-dup,]
			}
		

		if(nrow(x) >= 2) {
			o <- apply(combn(1:nrow(x),2),2,function(y){
				
				z <- x[y,]
				location <- c(z$qend[1]+1,z$qstart[2]-1);
				size <- diff(location) + 1
				#primers are in +for <-> -rev or -for <-> +rev order
				p <-  size <= size.range[2] && size >= size.range[1] && (all(unlist(z[,c("ptype","dir")]) == c(1,-1,1,-1)) | all(unlist(z[,c("ptype","dir")]) == c(-1,1,1,-1)))
				if(p) {
					
					list(
						dir=z$ptype[1],
						blastres=z,
						location=location,
						size=size
						)
					} else {
					list()
					}
				
				})
			o[unlist(lapply(o,length)) > 0]
			
			}
		
		
		
		})
	structure(o[unlist(lapply(o,length)) > 0],query=attr(x,"query"),class="insilicoPCR")
	
	
	}



insilico.tab <- function(x,blast.fields=NULL){
	ldply(lapply(x,function(x){
		ldply(lapply(x,function(x){
			b <- x$blastres;
			with(x,{
				a <- data.frame(size,dir,t(as.matrix(location)))
				names(a)[3:4] <- c("start","end")
				if(!is.null(blast.fields)) cbind.data.frame(a,t(as.matrix(unlist(b[,blast.fields])))) else a
				
				})
			}))
		}),.id="seq")
	}


#' Extract amplicons from sequences from \code{insilicoPCR} object
#' @param x object of class \code{insilicoPCR}
#' @return object of class \code{DNAStringSet}
#' @examples
#' ##Run local primer blast
#'
#' query <- readDNAStringSet("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/010/525/GCA_000010525.1_ASM1052v1/GCA_000010525.1_ASM1052v1_genomic.fna.gz")
#'
#' primers <- c("GGCTGGATCACCTCCTT","GCCWAGGCATCCDCC")
#' db <- makeprimerdb(primers)
#'
#' #Define BLAST arguments: these are the same that can be found in `RinsilicoPCR::blast.params$primers`
#' args <- list(
#' 		outfmt="'6 qacc sacc pident length qlen slen qstart sstart qend send evalue'",
#' 		max_target_seqs=20,
#' 		word_size=8
#' 		)
#'
#' blastres <- blast(query,db,args)
#' 
#' pairs <- primer.insilico(blastres)
#' 
#' amplicons <- extract(pairs)
#'
#' @export
extract <- function(x){
	if(!length(x)) return(NULL)
	fn <- attr(x,"query");
	dna <- readDNAStringSet(fn);
	names(dna) <- gsub(" .*","",names(dna));
	t <- insilico.tab(x);
	loc <- DNAStringSet(apply(t[,c("seq","start","end","dir")],1,function(x){
		s <- dna[[x['seq']]][x['start']:x['end']];
		if(x['dir'] == -1) reverseComplement(s) else s
		}))
		
	names(loc) <- with(t,{paste(seq,start,end,sep=".")})
	loc
	}	

