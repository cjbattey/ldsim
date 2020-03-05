library(ggplot2);library(reshape);library(cowplot)
theme_set(theme_classic()+theme(axis.line = element_blank(),
                                axis.text.y=element_blank(),
                                axis.ticks.y=element_blank(),
                                axis.title.y=element_blank(),
                                strip.background = element_blank()))

haplotype_freqs <- function(pd) {
  hf <- summary(factor(paste(pd$l1,pd$l2)))/sum(summary(factor(paste(pd$l1,pd$l2))))
  hf <- data.frame(haplotype=names(hf),frequency=hf)
  
  rownames(hf) <- NULL
  return(hf)
}

allele_freqs <- function(pd){
  af1 <- summary(factor(pd$l1))/sum(summary(factor(pd$l1)))
  af2 <- summary(factor(pd$l2))/sum(summary(factor(pd$l2)))
  af <- c(af1,af2)
  af <- data.frame(allele=names(af),frequency=af,locus=c(1,1,2,2))
  rownames(af) <- NULL
  return(af)
}

getLD <- function(allele1,allele2,pop,r2=FALSE){
  af1 <- sum(pop$l1==allele1)/nrow(pop)
  af2 <- sum(pop$l2==allele2)/nrow(pop)
  gf <- sum(pop$l1==allele1 & pop$l2==allele2)/nrow(pop)
  D <- gf-af1*af2
  if(r2){
    r2 <- D^2/(af1*(1-af1)*af2*(1-af2))
    return(r2)
  } else {
    return(D)
  }
}

plot_population <- function(pop,r){
  df <- data.frame(x=-0.02,xend=0.52,y=1:nrow(pop),yend=1:nrow(pop))
  pd <- cbind(df,pop)
  hf <- haplotype_freqs(pop)
  mh <- c("a b","a B","A b","A B")[!c("a b","a B","A b","A B") %in% hf$haplotype]
  if(length(mh)>0) hf <- rbind(hf,data.frame(haplotype=mh,frequency=0))
  
  af <- allele_freqs(pop)
  pd$ind <- unlist(lapply(1:(nrow(pd)/2),function(e) rep(e,2)))
  inds <- ddply(pd,.(ind),function(e) e[1,])
  inds$yend <- inds$y+1
  inds$x <- inds$x-0.01
  
  if(any(pd$recomb==1)){
    cplot <- ggplot(data=pd,aes(x=x,xend=xend,y=y,yend=yend))+
      ggtitle("Population")+
      xlab("recombination distance (r)")+
      ylab("Chromosomes")+
      theme(axis.title.y=element_text(angle=90,size=13),
            axis.text.x=element_text(size=12),
            axis.title.x=element_text(size=13))+
      geom_segment()+
      geom_segment(data=inds,aes(x=x,xend=x,y=y,yend=yend),color="grey",lwd=2)+
      scale_fill_manual(values = c("orange","red3","cornflowerblue","darkblue"),guide=F)+
      geom_segment(data=subset(pd,recomb==1)[seq(1,nrow(pd),2),],aes(x=r/2-r*0.25,xend=r/2+r*0.25,y=y,yend=y),color="white",lwd=2)+
      geom_segment(data=subset(pd,recomb==1)[seq(2,nrow(pd),2),],aes(x=r/2-r*0.25,xend=r/2+r*0.25,y=y,yend=y),color="white",lwd=2)+
      geom_segment(data=subset(pd,recomb==1)[seq(1,nrow(pd),2),],aes(x=r/2-r*0.25,xend=r/2+r*0.25,y=y,yend=y+1),color="purple",linetype=2)+
      geom_segment(data=subset(pd,recomb==1)[seq(2,nrow(pd),2),],aes(x=r/2-r*0.25,xend=r/2+r*0.25,y=y,yend=y-1),color="purple",linetype=2)+
      geom_label(aes(x=0,y=y,label=l1,fill=l1),color="white")+
      geom_label(aes(x=r,y=y,label=l2,fill=l2),color="white")
  } else {
    cplot <- ggplot(data=pd,aes(x=x,xend=xend,y=y,yend=yend))+
      ggtitle("Population")+
      xlab("recombination distance (r)")+
      ylab("Chromosomes")+
      theme(axis.title.y=element_text(angle=90),
            axis.text.x=element_text(size=12),
            axis.title.x=element_text(size=13))+
      geom_segment()+
      geom_segment(data=inds,aes(x=x,xend=x,y=y,yend=yend),color="grey",lwd=2)+
      scale_fill_manual(values = c("orange","red3","cornflowerblue","darkblue"),guide=F)+
      #geom_segment(data=subset(pd,recomb==1),aes(x=r/2,xend=r/2,y=y-0.25,yend=y+0.25),color="purple")+
      geom_label(aes(x=0,y=y,label=l1,fill=l1),color="white")+
      geom_label(aes(x=r,y=y,label=l2,fill=l2),color="white")
  }
  
  hfplot <- ggplot(data=hf,aes(x=haplotype,y=frequency))+
    theme_classic()+theme(axis.text=element_text(size=12),
                          axis.title=element_text(size=13))+
    ylim(0,1)+
    ggtitle("Haplotype Frequencies")+
    geom_bar(stat="identity")

  afplot <- ggplot(data=af,aes(x=allele,y=frequency))+
    facet_wrap(~locus,scales = "free_x",labeller = function(e) lapply(e,function(e) paste("locus",e)))+
    theme_classic()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=13))+
    ylim(0,1)+
    ggtitle("Allele Frequencies")+
    scale_fill_manual(values = c("orange","red3","cornflowerblue","darkblue"),guide=F)+
    geom_bar(stat="identity",aes(fill=allele))

  ggdraw()+
    draw_plot(cplot,0,0,0.6,1)+
    draw_plot(hfplot,0.6,0.5,0.4,0.5)+
    draw_plot(afplot,0.6,0,0.4,0.5)
}

makenewpop <- function(n){
  pop <- data.frame(l1=sample(c("A","a"),n,replace = T),
                    l2=sample(c("B","b"),n,replace = T))
  pop$recomb <- 0
  return(pop)
}

makegametes <- function(pop,r){ 
  newgen <- data.frame(l1=NA,l2=NA,recomb=0)[0,]
  chrs <- 1:nrow(pop)
  for(i in 1:(nrow(pop)/2)){
    ind <- sample(chrs,2)
    chrs <- chrs[!chrs %in% ind]
    parent <- pop[ind,]
    recomb <- rbinom(1,1,r)
    parent$recomb <- recomb
    if(recomb){
      tmp <- parent$l1[1]
      parent$l1[1] <- parent$l1[2]
      parent$l1[2] <- tmp
    }
    newgen <- rbind(newgen,parent)
  }
  return(newgen)
}

LDpred <- function(d1,r,gen){
  dpred <- c(d1)
  for(i in 2:(gen)){
    dpred[i] <- dpred[i-1]-dpred[i-1]*r
  }
  return(dpred)
}
