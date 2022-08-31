library(dplyr)
library(ggplot2)

args = commandArgs(TRUE)
if(length(args)==0){
  message('Arguments:')
  message("\t1: tsv with paths to identity and csv files for each sample. Columns should be 'sample', 'identity', 'length'")
  message('\t2: output PDF')
  message('\t3: output tsv')
}

## read paths
files.df = read.table(args[1], sep='\t', header=TRUE, as.is=TRUE)

## read identity csv files
id.df = lapply(1:nrow(files.df), function(ii){
  df = read.csv(files.df$identity[ii], as.is=TRUE, header=FALSE)
  colnames(df) = c('identity', 'n')
  df %>% mutate(identity=round(identity, 4)) %>%
    group_by(identity) %>% summarize(n=sum(n)) %>%
    mutate(sample=files.df$sample[ii])
}) %>% bind_rows

## summary table
id.sum = id.df %>%
  group_by(sample) %>% 
  arrange(identity) %>%
  mutate(prop=n/sum(n),
         prop.c=cumsum(prop),
         identity.geq.95=identity>=.95) %>% 
  summarize(mean.identity=sum(n*identity)/sum(n),
            median.identity=identity[which.min(abs(prop.c-.5))],
            prop.identity.geq.95=sum(n[identity.geq.95])/sum(n))

## width of the bar or point
id.width = .001
id.breaks = seq(0,1,id.width)
id.labels = id.breaks[-1] - id.width/2
id.b = id.df %>% mutate(identity.c=cut(identity, id.breaks, id.labels, include.lowest=TRUE),
                        identity.c=as.numeric(as.character(identity.c))) %>% 
  group_by(sample, identity.c) %>% summarize(n=sum(n)) %>% 
  group_by(sample) %>% arrange(identity.c) %>% mutate(prop=n/sum(n), prop.c=cumsum(prop))

id.sum = id.b %>% group_by(sample) %>%
  summarize(peak.identity=identity.c[which.max(prop)]) %>%
  merge(id.sum)

ggp.line.samp = ggplot(id.b, aes(x=identity.c, y=prop)) + 
    geom_area(alpha=.7) + 
    theme_bw() + 
    scale_color_brewer(palette="Set1") + 
    facet_grid(sample~.) + 
    theme(axis.text.y=element_blank(),
         strip.text.y=element_text(angle=0)) + 
    xlab('identity') + ylab('proportion of reads')

ggp.line.samp.70 = id.b %>% filter(identity.c>=.7) %>% 
    ggplot(aes(x=identity.c, y=prop)) + 
    geom_area(alpha=.7) + 
    theme_bw() + 
    scale_color_brewer(palette="Set1") + 
    facet_grid(sample~.) + 
    theme(axis.text.y=element_blank(),
         strip.text.y=element_text(angle=0)) + 
    xlab('identity') + ylab('proportion of reads')

## violin plot using sampled distribution
sampd <- function(df, N, coln='identity'){
    df = tibble(sample=df$sample[1], x=sample(df[,coln, drop=TRUE], N, replace=TRUE, prob=df$n/sum(df$n)))
    colnames(df)[2] = coln
    df
}
winsor <- function(x, u=NULL, l=NULL){
    if(any(x>u)) x[x>u] = u
    if(any(x<l)) x[x<l] = l
    x
}
id.s = id.df %>% group_by(sample) %>% do(sampd(., 10000))

ggp.id.violin = id.s %>% 
  ggplot(aes(x=sample)) +
  geom_violin(aes(y=winsor(identity, l=.7)), alpha=.6, fill='grey90') + 
  scale_fill_brewer(palette='Set1') + 
  geom_point(aes(y=median.identity), data=id.sum) + 
  xlab('sample') + ylab('identity') + 
  scale_y_continuous(breaks=seq(0,1,.1)) + 
  theme_bw() + 
  theme(legend.title=element_blank(), legend.position='top',
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

## read length csv files
len.df = lapply(1:nrow(files.df), function(ii){
  df = read.csv(files.df$length[ii], as.is=TRUE, header=FALSE)
  colnames(df) = c('length', 'n')
  df %>% mutate(sample=files.df$sample[ii])
}) %>% bind_rows

len.df = len.df %>% group_by(sample) %>% 
  arrange(length) %>%
  mutate(prop=n/sum(n),
         prop.c=cumsum(prop),
         length.log10=log10(length)) %>%
  arrange(desc(length)) %>%
  mutate(nxx=cumsum(n*length/1e6)/sum(n*length/1e6))


ggp.nxx = len.df %>% 
    group_by(sample) %>% sample_n(1000) %>% 
    ungroup %>% arrange(desc(sample)) %>%
    ggplot(aes(x=nxx*100, y=length/1e3, colour=sample)) +
    geom_vline(xintercept=50, linetype=2) +
    scale_x_continuous(breaks=seq(0,100,10)) + 
    geom_line(alpha=7) + theme_bw() +
    xlab('cumulative sequence percent') +
    ylab('read length (Kbp)') + 
    scale_color_brewer(palette='Set1') + 
    theme(legend.position=c(.99, .99), legend.justification = c(1,1),
         legend.title=element_blank())

ggp.nxx.95 = len.df %>% filter(nxx>=.05) %>% 
    group_by(sample) %>% sample_n(1000) %>% 
    ungroup %>% arrange(desc(sample)) %>%
    ggplot(aes(x=nxx*100, y=length/1e3, colour=sample)) +
    geom_vline(xintercept=50, linetype=2) +
    scale_x_continuous(breaks=seq(0,100,10)) + 
    geom_line(alpha=7) + theme_bw() +
    xlab('cumulative sequence percent') +
    ylab('read length (Kbp)') + 
    scale_color_brewer(palette='Set1') + 
    theme(legend.position=c(.99, .99), legend.justification = c(1,1),
         legend.title=element_blank())

## log-scaled distribution
len.width = (1+max(len.df$length.log10)) / 1000
len.breaks = seq(0,max(len.df$length.log10)+1,len.width)
len.labels = len.breaks[-1] - len.width/2
len.b = len.df %>% mutate(length.c=cut(length.log10, len.breaks, len.labels, include.lowest=TRUE),
                          length.c=as.numeric(as.character(length.c))) %>% 
  group_by(sample, length.c) %>% summarize(n=sum(n)) %>% 
  group_by(sample) %>% arrange(length.c) %>% mutate(prop=n/sum(n), prop.c=cumsum(prop))

ggp.len.dist.log10 = ggplot(len.b, aes(x=length.c, y=prop, colour=sample)) + 
    geom_line(size=1) + 
    theme_bw() + 
    scale_color_brewer(palette='Set1') + 
    xlab('log10(read length)') + ylab('proportion of reads')

## summary table
len.sum = len.df %>%
  mutate(length.geq.10kb=length>=1e4) %>%
  group_by(sample) %>% 
  summarize(mean.length=sum(n*length)/sum(n),
            median.length=length[which.min(abs(prop.c-.5))],
            n50=length[which.min(abs(nxx-.5))],
            total.gbp=sum(n*length)/1e9,
            prop.length.geq.95=sum(n[length.geq.10kb])/sum(n))

## save graphs in pdf
pdf(args[2], 9, 5)
print(ggp.line.samp)
print(ggp.line.samp.70)
print(ggp.id.violin)
print(ggp.len.dist.log10)
print(ggp.nxx)
print(ggp.nxx.95)
dev.off()

## save sum in tsv
write.table(merge(id.sum, len.sum), file=args[3], sep=',', row.names=FALSE, quote=FALSE)
