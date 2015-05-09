
### copied and modified on 10 Apr, 2015 to generate testing cases for discrete
### dynamic bayesian network, mainly as model of gene regulatory
### network with time delays.

parent.index <- function(states, p.factors) {
  # states are the states of the parents (0-based small integers),
  # p.factors are the corresponding factors,
  # calculate an index to represent the vector of states
  sum(states*p.factors)
}

parent.states <- function(idx, p.states) {
  # basically inverse of parent.index, to turn index back into vector
  # of states, but need the number of states of the parents
  if(length(p.states) <= 0)
    0
  else {
    r <- rep(0, length(p.states));
    for(i in 1:length(p.states)) {
      r[i] <- idx %% p.states[i];
      idx <- idx %/% p.states[i];
    }
    r
  }
}

gen.dist <- function(states, p.bias) {
  ## states is the number of states for the current node,
  ## the most probable state has probability p.bias.
  ##  r <- runif(states);
  r <- rep((1-p.bias)/(states-1), states);
  s <- sample(1:states, 1);
  r[s] <- p.bias;
  r
}

gen.parents.dist <- function(ps.idx, ps.delays, p.states, states, p.bias) {
  ## get the CPT and other things, once the parents (ps.idx) and delays (ps.delays) are decided.
  ## states is the number of states for the current node,
  ## p.states are number of states for the parents in ps.idx
  ## p.bias is the probability of the most probable state in each conditional distribution.

  ## prepare the multiplcative factors
  p.f <- cumprod(c(1, head(p.states, -1)));
  ## conditional distribution table is a long vector, where every
  ## states numbers represent a conditional distribution, and the
  ## states of the parents are compacted into an index.
  n.p.states <- prod(p.states);
  p.cpt <- rep(0, states*n.p.states);
  j <- 0;
  for(i in 1:n.p.states) {
    p.cpt[j + 1:states] <- gen.dist(states, p.bias);
    j <- j+states;
  }
  #
  list(states=states, parents=ps.idx, delays=ps.delays, p.states=p.states, p.factors=p.f, n.p.states=n.p.states, p.cpt=p.cpt)
}

gen.parents <- function(states, n, n.states, n.parents, max.delay, p.bias) {
  ## randomly select parents, and generate the CPT
  ## states is the number of states for the current node,
  ## n.states are number of states for all n nodes.
  ## p.bias is the probability of the most probable state in each conditional distribution.
  ps <- sample(1:n, n.parents); # indices of parents
  ds <- sample(1:max.delay, n.parents, replace=TRUE); # delays for each parent
  p.states <- n.states[ps]; # number of states for the parents

  gen.parents.dist(ps, ds, p.states, states, p.bias)
}

gen.small.hidden <- function(m,n, n.states, max.delay, p.bias) {
  ## m nodes as parent, of a hidden node, then n children of the hidden node
  ## the hidden node is placed as the last
  ## n.states is the number of states of each node.
  ng <- m+n+1;
  r <- list();

  ds <- sample(1:max.delay, m+n, replace=TRUE);
  ## m parent nodes of the hidden node
  if(m > 0) {
    for(i in 1:m) {
      r[[i]] <- gen.parents.dist(integer(0), integer(0), integer(0), n.states, p.bias);
    }
  }
  ## n children nodes of the hidden node
  for(i in m+(1:n)) {
    r[[i]] <- gen.parents.dist(ng, ds[i], n.states, n.states, p.bias);
  }
  ## the hidden node
  if(m > 0) {
    r[[ng]] <- gen.parents.dist(1:m, ds[1:m], rep(n.states,m), n.states, p.bias);
  } else {
    r[[ng]] <- gen.parents.dist(NULL, NULL, rep(n.states,m), n.states, p.bias);
  }
  ##
  r
}

gen.grn <- function(n, n.states, possible.parents, max.delay, p.bias) {
  # n is the number of genes (nodes)
  # possible.parents is the vector of possible number of parents to choose from.
  # n.states is the number of states of each gene, so should be of length n. It will be repeated as necessary to length n.
  n.states <- rep(n.states, length.out=n);
  r <- list();
  mp <- sample(possible.parents, n, replace=TRUE);
  for(i in 1:n) {
    r[[i]] <- gen.parents(n.states[i], n, n.states, mp[i], max.delay, p.bias);
  }
  r
}

gen.exp <- function(grn, nps, max.delay) {
  # simulate the discrete expression of n genes in grn, with nps points.
  # returns a matrix of nps by n.
  # the states are non-negative small integers.
  n <- length(grn);
  r <- matrix(0, nrow=nps+max.delay, ncol=n);
  # generate the initial expression, uniformly
  for(i in 1:n) {
    r[1:max.delay,i] <- sample(0:(grn[[i]]$states - 1), max.delay, replace=TRUE);
  }
  #
  # print(r[1:max.delay,]);
  #
  for(i in (max.delay+1):(max.delay+nps)) {
    # cat("time ",i,"\n");
    for(j in 1:n) {
      # cat("gene",j,"\n");
      # get the states of parents
      ns <- grn[[j]]$states;
      ps <- grn[[j]]$parents;
      ds <- grn[[j]]$delays;
      if(length(ps) > 0) {
        s <- rep(0, length(ps));
        for(k in 1:length(ps)) {
          s[k] <- r[i-ds[k], ps[k]];
        }
        # generate the conditional state
        idx <- parent.index(s, grn[[j]]$p.factors);
        # cat(s, " idx:",idx,"\n");
        # cat("prob: ", grn[[j]]$p.cpt[idx*ns + (1:ns)], "\n");
        c.prob <- grn[[j]]$p.cpt[idx*ns + (1:ns)];
      } else {
        c.prob <- grn[[j]]$p.cpt[1:ns];
      }
      r[i,j] <- sample(0:(ns-1), 1, prob=c.prob);
    }
  }
  #
  r[(max.delay+1):(max.delay+nps),]
}

##
print.exp <- function(d) {
  n <- nrow(d);
  for(i in 1:n) {
    cat(d[i,], sep=" ");
    cat("\n");
  }
}

## output the GRN
print.grn.links <- function(grn) {
  # output links and delays
  # all indices 0-based
  n <- length(grn); # number of genes
  for(i in 1:n) {
    ps <- grn[[i]]$parents;
    ds <- grn[[i]]$delays;
    if(length(ps) > 0) {
      for(j in 1:length(ps)) {
        cat("To:",i-1, "From:",ps[j]-1, "Delay:",ds[j]);
        cat("\n");
      }
    }
  }
}

## output also the conditional distributions
print.grn <- function(grn) {
  # output in a raw form
  # all indices 1-based
  n <- length(grn); # number of genes
  for(i in 1:n) {
    cat("=== Gene", i, "\n");
    ns <- grn[[i]]$states;
    cat("state:", ns, "\n");
    cat("parents: "); cat(grn[[i]]$parents, sep=" "); cat("\n");
    cat("delays: "); cat(grn[[i]]$delays, sep=" "); cat("\n");
    p.s <- grn[[i]]$p.states;
    cat("parents states: "); cat(p.s, sep=" "); cat("\n");
    n.ps <- grn[[i]]$n.p.states;
    cat("CPT:\n");
    for(j in 1:n.ps) {
      # parent states config : conditional distribution
      cat(parent.states(j-1, p.s), sep=" ");
      cat(" : ");
      cat(grn[[i]]$p.cpt[(j-1)*ns + 1:ns], sep=" ");
      cat("\n");
    }
    cat("=== End Gene", i, "\n");
  }
}

####
gen.one.small.case <- function(out.prefix, m,n, nps, segs, n.states=3, max.delay=4, p.bias=0.6) {
  ng <- m + n; ## number of observed genes, the last one is hidden
  
  grn <- gen.small.hidden(m,n, n.states, max.delay, p.bias);

  # grn in simple format
  sink(paste(out.prefix, "_grn.txt", sep=""));
  print.grn.links(grn);
  sink();

  # grn with CPT
  sink(paste(out.prefix, "_cpt.txt", sep=""));
  print.grn(grn);
  sink();

  # expression data, one long segment
  d <- gen.exp(grn, max(nps), max.delay);
  for(np in nps) {
    ## hidden
    sink(paste(out.prefix, "_nps",np,".txt", sep=""));
    print.exp(d[1:np,1:ng]);
    sink();
    ## complete
    sink(paste(out.prefix, "_nps",np,"_complete.txt", sep=""));
    print.exp(d[1:np,]);
    sink();
  }

  # multiple short segments
  for(i in 1:length(segs)) {
    d <- gen.exp(grn, segs[i], max.delay);
    ## hidden
    sink(paste(out.prefix, "_s",i,".txt", sep=""));
    print.exp(d[,1:ng]);
    sink();
    ## complete
    sink(paste(out.prefix, "_s",i,"_complete.txt", sep=""));
    print.exp(d);
    sink();
  }
}

gen.small.hidden.cases <- function() {
  out <- "small";
  md <- 4;
  ## need only generate the large number that we will use, then take
  ## fragments if we want to test shorter series
  nps <- c(25, 50, 100, 200, 400, 800);
  segs <- c(29, 31, 31, 31, 32,  28, 32, 32, 30, 29,
            23, 28, 28, 25, 28,  26, 27, 30, 33, 27,
            28, 25, 29, 29, 31,  29, 30, 24, 29, 26,
            27, 31);
  dir.create(out, mode="0755");
  for(m in c(0,1,2,3)) { ## number of parents
    for(n in c(2,3,4,5)) { ## number of children
      outp <- paste(out,"/m",m,"n",n, sep="");
      dir.create(outp, mode="0755");
      for(b in c(0.65,0.75,0.85)) { ## bias
        for(r in 1:20) {
          name <- paste(outp,"/b",b,"r",r, sep="");
          gen.one.small.case(name, m,n,nps, segs, max.delay=md, p.bias=b);
          cat(name,"\n");
        }
      }
    }
  }
}

#######################################################################################
read.one.grn <- function(ng, name) {
  ## since the GRN may not have connection for all genes, so the number of genes is given.
  ## read GRN in this format: (0-based)
  ## To: 0 From: 2 Delay: 1
  ## To: 1 From: 2 Delay: 1
  ## To: 2 From: 1 Delay: 4
  tmp <- try(read.table(name));
  rD <- matrix(0, nrow=ng,ncol=ng);
  if(! inherits(tmp, 'try-error')) {
    for(i in 1:nrow(tmp)) {
      from <- 1 + tmp[i,4];
      to <- 1 + tmp[i,2];
      if(tmp[i,6] > 0) {
        rD[from,to] <- tmp[i,6];
      }
    }
  }
  rD
}

delays.to.grn <- function(delays, n.states, p.bias) {
  ## delays[from,to] = delay;
  ## n.states is the number of states of each node.
  ng <- nrow(delays);

  r <- list();
  for(i in 1:ng) {
    idx <- which(delays[,i] > 0);
    ds <- delays[idx,i];
    r[[i]] <- gen.parents.dist(idx, ds, rep(n.states, length(idx)), n.states, p.bias);
  }
  ##
  r
}

gen.one.small.nh.case <- function(out.prefix, tmp.grn, m,n, nps, segs, n.states=3, max.delay=4, p.bias=0.6) {
  ng <- m + n; ## number of observed genes

  ## try not to have empty GRN
  if(sum(tmp.grn > 0) == 0) {
    grn <- gen.small.hidden(m,n-1, n.states, max.delay, p.bias);
  } else {
    grn <- delays.to.grn(tmp.grn,  n.states, p.bias);
  }

  # grn in simple format
  sink(paste(out.prefix, "_grn.txt", sep=""));
  print.grn.links(grn);
  sink();

  # grn with CPT
  sink(paste(out.prefix, "_cpt.txt", sep=""));
  print.grn(grn);
  sink();

  # expression data, one long segment
  # all are complete
  d <- gen.exp(grn, max(nps), max.delay);
  for(np in nps) {
    ## complete
    sink(paste(out.prefix, "_nps",np,".txt", sep=""));
    print.exp(d[1:np,]);
    sink();
  }

  # multiple short segments
  for(i in 1:length(segs)) {
    d <- gen.exp(grn, segs[i], max.delay);
    ## complete
    sink(paste(out.prefix, "_s",i,".txt", sep=""));
    print.exp(d);
    sink();
  }
}

gen.small.non.hidden.cases <- function() {
  out <- "small_nh";
  md <- 4;
  ## need only generate the large number that we will use, then take
  ## fragments if we want to test shorter series
  nps <- c(25, 50, 100, 200, 400, 800);
  segs <- c(29, 31, 31, 31, 32,  28, 32, 32, 30, 29,
            23, 28, 28, 25, 28,  26, 27, 30, 33, 27,
            28, 25, 29, 29, 31,  29, 30, 24, 29, 26,
            27, 31);
  dir.create(out, mode="0755");
  for(m in c(0,1,2,3)) { ## number of parents
    for(n in c(2,3,4,5)) { ## number of children
      outp <- paste(out,"/m",m,"n",n, sep="");
      dir.create(outp, mode="0755");
      for(b in c(0.65,0.75,0.85)) { ## bias
        for(r in 1:20) {
          name <- paste(outp,"/b",b,"r",r, sep="");
          tmp.grn <- read.one.grn(m+n, paste("small/mit_fast_h_grn/m",m,"n",n,"/b",b,"r",r,"_nps800_out_grn.txt", sep=""));
          gen.one.small.nh.case(name, tmp.grn, m,n,nps, segs, max.delay=md, p.bias=b);
          cat(name,"\n");
        }
      }
    }
  }
}

#######################################################################################
random.fill <- function(L, es) {
  ## fill in es for length L as much as possible, the rest are randomly picked
  if(L < length(es)) {
    sample(es,L, replace=FALSE)
  } else {
    r <- rep(es, L %/% length(es));
    n <- L %% length(es);
    if(n > 0) {
      r <- c(r, sample(es, n, replace=FALSE));
    }
    sample(r,L, replace=FALSE)
  }
}

gen.misc.hidden <- function(ng, nh, max.parents, max.delay) {
  # ng observed nodes (the first ng), nh hidden nodes.
  # ng should be large enough.
  # each hidden node will have 0,1,2 or 3 parents, and 2,3,4 or 5 children, all distinct.
  # the hidden node is placed as the last
  # Returns a matrix of (ng+nh) by (ng+nh), for the delays[from,to].
  n <- ng + nh;
  rdelay <- matrix(0, nrow=n, ncol=n);
  ps <- c(0,1,2,3);
  cs <- c(2,3,4,5);
  np <- random.fill(nh, ps); ## number of parents for each hidden node
  nc <- random.fill(nh, cs); ## number of children for each hidden node
  nps <- sum(np);
  ncs <- sum(nc);
  n.other <- ng - nps - ncs;
  mp <- sample(1:max.parents, ng - ncs, replace=TRUE); ## number of parents
  for(i in 1:(ng - ncs)) { ## fill in the columns one by one, (i,j) is i --> j
    if(i <= n.other) { ## others, can use any observed gene as parents
      idx <- sample(1:ng, mp[i]);
    } else { ## parents of hidden nodes, parents cannot be children of hidden nodes
      idx <- sample(1:(ng-ncs), mp[i]);
    }
    rdelay[idx,i] <- sample(1:max.delay, mp[i], replace=TRUE);
  }
  ## the hidden nodes
  pj <- n.other;
  cj <- ng - ncs;
  for(i in 1:nh) {
    ## parents
    if(np[i] > 0) {
      idx <- pj + (1:np[i]);
      rdelay[idx,ng+i] <- sample(1:max.delay, np[i], replace=TRUE);
    }
    ## children, all nc > 0
    idx <- cj + (1:nc[i]);
    rdelay[ng+i,idx] <- sample(1:max.delay, nc[i], replace=TRUE);
    ##
    pj <- pj + np[i];
    cj <- cj + nc[i];
  }
  ##
  rdelay
}

gen.one.misc.case <- function(out.prefix, ng,nh, nps, segs, n.states=3, max.parents=2, max.delay=4, p.bias=0.6) {
  tmp <- gen.misc.hidden(ng,nh, max.parents, max.delay);
  grn <- delays.to.grn(tmp, n.states, p.bias);

  # grn in simple format
  sink(paste(out.prefix, "_grn.txt", sep=""));
  print.grn.links(grn);
  sink();

  # grn with CPT
  sink(paste(out.prefix, "_cpt.txt", sep=""));
  print.grn(grn);
  sink();

  # expression data, one long segment
  d <- gen.exp(grn, max(nps), max.delay);
  for(np in nps) {
    ## hidden
    sink(paste(out.prefix, "_nps",np,".txt", sep=""));
    print.exp(d[1:np,1:ng]);
    sink();
    ## complete
    sink(paste(out.prefix, "_nps",np,"_c.txt", sep=""));
    print.exp(d[1:np,]);
    sink();
  }

  # multiple short segments
  for(i in 1:length(segs)) {
    d <- gen.exp(grn, segs[i], max.delay);
    ## hidden
    sink(paste(out.prefix, "_s",i,".txt", sep=""));
    print.exp(d[,1:ng]);
    sink();
    ## complete
    sink(paste(out.prefix, "_s",i,"_c.txt", sep=""));
    print.exp(d);
    sink();
  }
}

gen.misc.cases <- function() {
  out <- "misc";
  md <- 4;
  ## need only generate the large number that we will use, then take
  ## fragments if we want to test shorter series
  nps <- c(25, 50, 100, 200, 400, 800, 1000, 1200, 1400, 1600);
  segs <- c(29, 31, 31, 31, 32,  28, 32, 32, 30, 29,
            23, 28, 28, 25, 28,  26, 27, 30, 33, 27,
            28, 25, 29, 29, 31,  29, 30, 24, 29, 26,
            27, 31,
            29, 31, 31, 31, 32,  28, 32, 32, 30, 29,
            23, 28, 28, 25, 28,  26, 27, 30, 33, 27,
            28, 25, 29, 29, 31,  29, 30, 24, 29, 26,
            27, 31);
  dir.create(out, mode="0755");
  for(ng in c(50, 100)) { # number of observed genes
    ## with 10% hidden nodes
    nh <- ng/10;
    outp <- paste(out,"/n",ng, sep="");
    dir.create(outp, mode="0755");
    for(b in c(0.65,0.75,0.85)) { ## bias
      dir.create(paste(outp,"/b",b, sep=""), mode="0755");
      for(r in 1:40) {
        name <- paste(outp,"/b",b,"/r",r, sep="");
        gen.one.misc.case(name, ng,nh,nps, segs, max.delay=md, p.bias=b);
        cat(name,"\n");
      }
    }
  }
}

####
