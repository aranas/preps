function div = finddivisor(a)

K=1:a;
div = K(rem(a,K)==0);