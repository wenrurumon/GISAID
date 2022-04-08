
library(phytools)
library(ggplot2)
setwd('/Users/wenrurumon/Documents/postdoc/covid19/fasta')

raw <- readLines('gisaid_variants_statistics.tsv')
raw <- do.call(rbind,strsplit(raw,'\t'))
colnames(raw) <- raw[1,]
raw <- as.data.table(raw[-1,])
raw <- raw %>% 
  filter(Type=='Lineage') %>% 
  group_by(date=`Week prior to`,lineage=Value) %>% summarise(n=sum(as.numeric(`Submission Count`))) %>%
  mutate(date=as.Date(date))
raw$group <- out$group[match(raw$lineage,out$lineage)]
raw <- raw %>%mutate(group=ifelse(is.na(group),'Others',group))

temp <- raw %>% group_by(date,group) %>% summarise(n=sum(n))
temp <- temp %>% merge(temp %>% group_by(date) %>% summarise(p=sum(n))) %>% mutate(p=n/p) %>%
  filter(date>='2021-10-01')

x <- temp %>% acast(date~group,value.var='p',fill=0)
x <- apply(x,2,function(x){
  splinefun((1:length(x)-1)*7+1,x)(1:183)
})
rownames(x) <- paste(min(temp$date)+0:182)
x[x<0] <- 0
x <- x/rowSums(x)
x <- melt(x) %>% select(date=1,group=2,p=3) %>% mutate(date=as.Date(date))
ggplot() + geom_area(
  data=x,
  aes(x=date,y=p*100,fill=group)
) + scale_fill_manual(values=c(BA.1='#F50E01',BA.1.1='#41BF00',BA.2='#AD08E3',BA.3='#F78001',Delta='blue',Others='gray')) + 
  labs(x='',y='Prevalence (%)',fill='') + 
  scale_x_date(date_breaks='1 month',labels=date_format("%Y %b")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(size=15,face='bold'),
        legend.position='right')
ggsave("figb1.eps", width = 8.0, height = 5.0,device=cairo_ps)

tree <- read.tree("nramp.nwk")
t2drop <- c('BA.3_187','BA.3_151')
tree <- drop.tip(tree,t2drop)
tree$edge.length <- tree$edge.length*nchar(temp[2])
PatristicDistMatrix<-cophenetic.phylo(tree)
group_info <- data.frame(label=tree$tip.label,group=sapply(strsplit(tree$tip.label,'_'),function(x){x[1]}))
tree1 <- full_join(tree,group_info,by='label')
ggtree(tree1,aes(color=group),size=1) +
  geom_tippoint(size=2,aes(fill=group)) +
  geom_tiplab(size=0) +
  guides(color=guide_legend(override.aes=list(size=10))) +
  scale_colour_manual(values=c(BA.1='#F50E01',BA.1.1='#41BF00',BA.2='#AD08E3',BA.3='#F78001')) +
  theme_tree2(legend.position='none',text=element_text(size=15,face='bold'),
              panel.grid.major   = element_line(color='black', size=.2),
              panel.grid.minor   = element_line(color='grey', size=.2),
              panel.grid.major.y=element_blank(),
              panel.grid.minor.y=element_blank(),
              axis.line.x = element_line(size=0))
ggsave("figb2.eps", width = 8.0, height = 5.0,device=cairo_ps)
