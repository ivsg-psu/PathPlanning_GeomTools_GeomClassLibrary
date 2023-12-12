for ith_frame = 1:100
    drawframe(ith_frame);
    pause(0.01);
end

function drawframe(f)
rng default;
P = rand([round(f*.5)+2 2]);
DT = delaunayTriangulation(P)
IC = incenter(DT);
triplot(DT)
tri1 = DT.ConnectivityList(1, [1:end 1]);
hold on
plot(P(tri1,1), P(tri1,2), 'r') % first triangle is red
plot(IC(:,1),IC(:,2),'.c') % incenter dots
axis('square','equal','off');
hold off
end