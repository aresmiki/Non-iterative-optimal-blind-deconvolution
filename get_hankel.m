function data=get_hankel(rbt1,np)
    N=length(rbt1);
    data=[];
    for i=1:floor(N-np+1)
        data(i,:)=rbt1(i:(i+np-1));
    end
end