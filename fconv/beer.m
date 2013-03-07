function apod=beer(N,MD);

beer=zeros(N,1);
beer(1)=1;
for i=2:N/2
        if i <= MD
                beer(i)=(1-((i-1)/MD)^2)^2;
        else
                beer(i)=0;
        end
end
beer(N/2+1:N)=flipud(beer(1:N/2));

apod=beer;
