---
title: "Little Problems of Perl"
date: '2019-02-01'
slug: perlproblems
categories: ["data at fingertips"]
tags: ["Bioinformatics"]
---

### undefined symbol: Perl\_Gthr\_key\_ptr

How to solve the problem of "symbol lookup error: /home/zhengh42/perl5/lib/perl5/x86_64-linux-thread-multi/auto/Clone/Clone.so: undefined symbol: Perl\_Gthr\_key\_ptr"?

My specific solution:

Look at PERL5LIB environment variable: `perl -V` Or in ~/.bash_rc

Then reset the PERL5LIB variable: `unset PERL5LIB`

### Install perl packages locally
```
% cpan
cpan> o conf makepl_arg INSTALL_BASE=/mydir/perl
cpan> o conf commit
cpan> install themodule
```

