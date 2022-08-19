library(dplyr)
library(ggplot2)

args = commandArgs(TRUE)
if(length(args)==0){
  message('Arguments:')
  message('\t1: identity_distribution.csv')
  message('\t2: length_distribution.csv')
  message('\t3: output PDF')
  message('\t4: output tsv')
}

##
## identity
id.df = read.csv(args[1], as.is=TRUE, header=FALSE)
colnames(id.df) = c('identity', 'n')

## 4 digits precision is enough for the graphs
id.df = id.df %>% mutate(identity=round(identity, 4)) %>%
  group_by(identity) %>% summarize(n=sum(n))

## width of the bar or point
id.width = .001
id.breaks = seq(0,1,id.width)
id.labels = id.breaks[-1] - id.width/2
id.b = id.df %>% mutate(identity.c=cut(identity, id.breaks, id.labels, include.lowest=TRUE),
                        identity.c=as.numeric(as.character(identity.c))) %>% 
  group_by(identity.c) %>% summarize(n=sum(n)) %>% 
  ungroup %>% arrange(identity.c) %>% mutate(prop=n/sum(n), prop.c=cumsum(prop))

ggp.id.dist = ggplot(id.b, aes(x=identity.c, y=prop)) + 
    geom_line(size=1) + 
    theme_bw() + 
    xlab('identity') + ylab('proportion of reads')

## zoom to top 90% of data
ggp.id.dist.90 =
  id.b %>% filter(prop.c>.1) %>%
  ggplot(aes(x=identity.c, y=prop)) + 
  geom_line(size=1) + 
  theme_bw() + 
  xlab('identity') + ylab('proportion of reads')

## cumulative distribution curve
ggp.id.cdist = ggplot(id.b, aes(x=identity.c, y=prop.c)) + 
    geom_line(size=1) + 
    theme_bw() + 
    xlab('identity') + ylab('cumulative proportion of reads') + 
    scale_y_continuous(breaks=seq(0,1,.1))

## summary table
id.sum = id.df %>%
  arrange(identity) %>%
  mutate(prop=n/sum(n),
         prop.c=cumsum(prop),
         identity.geq.95=identity>=.95) %>% 
  summarize(mean.identity=sum(n*identity)/sum(n),
            median.identity=identity[which.min(abs(prop.c-.5))],
            prop.identity.geq.95=sum(n[identity.geq.95])/sum(n))

##
## length distribution
len.df = read.csv(args[2], as.is=TRUE, header=FALSE)
colnames(len.df) = c('length', 'n')
len.df = len.df %>%
  arrange(length) %>%
  mutate(prop=n/sum(n),
         prop.c=cumsum(prop),
         length.log10=log10(length)) %>%
  arrange(desc(length)) %>%
  mutate(nxx=cumsum(n*length)/sum(n*length))

## width of the bar or point
len.width = (1+max(len.df$length)) / 1000
len.breaks = seq(0,max(len.df$length)+1,len.width)
len.labels = len.breaks[-1] - len.width/2
len.b = len.df %>% mutate(length.c=cut(length, len.breaks, len.labels, include.lowest=TRUE),
                          length.c=as.numeric(as.character(length.c))) %>% 
  group_by(length.c) %>% summarize(n=sum(n)) %>% 
  ungroup %>% arrange(length.c) %>% mutate(prop=n/sum(n), prop.c=cumsum(prop))

## distribution
ggp.len.dist = ggplot(len.b, aes(x=length.c, y=prop)) + 
    geom_line(size=1) + 
    theme_bw() + 
    xlab('read length') + ylab('proportion of reads')

## cumulative distribution
ggp.len.cdist = ggplot(len.b, aes(x=length.c, y=prop.c)) + 
    geom_line(size=1) + 
    theme_bw() + 
    xlab('read length') + ylab('cumulative proportion of reads') + 
    scale_y_continuous(breaks=seq(0,1,.1))

## NXX graph
ggp.nxx = ggplot(len.df, aes(x=nxx, y=length)) +
  geom_vline(xintercept=.5, linetype=2) + 
  geom_line() + theme_bw() +
  xlab('cumulative sequence proportion') +
  ylab('read length')

## log-scaled distribution
len.width = (1+max(len.df$length.log10)) / 1000
len.breaks = seq(0,max(len.df$length.log10)+1,len.width)
len.labels = len.breaks[-1] - len.width/2
len.b = len.df %>% mutate(length.c=cut(length.log10, len.breaks, len.labels, include.lowest=TRUE),
                          length.c=as.numeric(as.character(length.c))) %>% 
  group_by(length.c) %>% summarize(n=sum(n)) %>% 
  ungroup %>% arrange(length.c) %>% mutate(prop=n/sum(n), prop.c=cumsum(prop))

ggp.len.dist.log10 = ggplot(len.b, aes(x=length.c, y=prop)) + 
    geom_line(size=1) + 
    theme_bw() + 
    xlab('log10(read length)') + ylab('proportion of reads')

## summary table
len.sum = len.df %>%
  mutate(length.geq.10kb=length>=1e4) %>% 
  summarize(mean.length=sum(n*length)/sum(n),
            median.length=length[which.min(abs(prop.c-.5))],
            n50=length[which.min(abs(nxx-.5))],
            total.gbp=sum(n*length)/1e9,
            prop.length.geq.95=sum(n[length.geq.10kb])/sum(n))

## save graphs in pdf
pdf(args[3], 9, 5)
print(ggp.id.dist)
print(ggp.id.dist.90)
print(ggp.id.cdist)
print(ggp.len.dist)
print(ggp.len.dist.log10)
print(ggp.len.cdist)
print(ggp.nxx)
dev.off()

## save sum in tsv
write.table(cbind(id.sum, len.sum), file=args[4], sep=',', row.names=FALSE, quote=FALSE)
