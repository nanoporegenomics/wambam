## TODO

## violin plot using sampled distribution
sampleDist <- function(N, lengths, ns){
  sample(lengths, N, replace=TRUE, prob=ns/sum(ns))
}

len.s = tibble(length=sampleDist(10000, len.df$length, len.df$n))
ggplot(len.s, aes(x=1, y=length)) +
  geom_violin() +
  theme_bw()
