while(1)
    for ith_frame = 1:25
    drawframe(ith_frame);
    pause(0.01);
    end
end
 
function drawframe(f)
if round((f-1)/4)==(f-1)/4 
D=findobj('Tag','SL')
for i=1:length(D)
    D(i).Visible='off';
end
ax=gca;
hold(ax,'on');
ax.XTick=[];
ax.YTick=[];
ax.XLim=[0,800];
ax.YLim=[0,600];
ax.DataAspectRatio=[1 1 1];
% =========================================================================
% 颜色预定义，注意此处是hsv格式
cClouds=[330,25,100];  % 云的颜色
cSky=[220,50,50];      % 天空的颜色
cFurther=[230,25,90];  % 远山的颜色
cCloser=[210,70,10];   % 近山的颜色
% =========================================================================
% 绘图函数调用
ax.Color=hsv2rgb(cFurther./[360,100,100]); % 背景为远山的颜色
drawSky(cSky,cFurther)                     % 画出天空颜色渐变效果 
drawClouds(cClouds)                        % 画出彩色云朵效果
drawMountains(cFurther,cCloser)            % 画出山脉效果
end
% =========================================================================
% 功能函数：
% -------------------------------------------------------------------------
% 渐变背景生成函数
    function drawSky(colSky,colFurther)
        % 颜色由hsv转rgb
        colSky=hsv2rgb(colSky./[360,100,100]);
        colFurther=hsv2rgb(colFurther./[360,100,100]);
        
        %构建渐变色网格
        [XMesh,YMesh]=meshgrid(1:800,301:600);
        ZMesh=zeros(size(XMesh));
        CMesh=vColorMat([800,300],[colFurther;colSky]);
        surf(XMesh,YMesh,ZMesh,'CData',CMesh,'EdgeColor','none','Tag','SL');
    end
% -------------------------------------------------------------------------
% 云绘制函数
    function drawClouds(colClouds)
        colClouds=hsv2rgb(colClouds./[360,100,100]);
        % 随机噪声生成
        [X,Y]=meshgrid(linspace(0,1,500));
        CLX=(-cos(X.*2.*pi)+1).^.2;
        CLY=(-cos(Y.*2.*pi)+1).^.2;
        r=(X-.5).^2+(Y-.5).^2;
        alp=abs(ifftn(exp(3i*rand(500))./r.^.8)).*(CLX.*CLY);
        alp=alp./max(alp,[],'all');
        
        CMesh=zeros([size(alp),3]);
        CMesh(:,:,1)=colClouds(1);
        CMesh(:,:,2)=colClouds(2);
        CMesh(:,:,3)=colClouds(3);
        
        % 越向下、云越透明
        dy=(1:500)./500.*0.8+0.2;
        image([0,800],[350,600],CMesh,'AlphaData',alp.*(dy'),'Tag','SL');
    end
% -------------------------------------------------------------------------
% 山峰绘制函数
    function drawMountains(colFurther,colCloser)
        [X,Y]=meshgrid(linspace(0,1,800));
        CLX=(-cos(X.*2.*pi)+1).^.2;
        CLY=(-cos(Y.*2.*pi)+1).^.2;
        r=(X-.5).^2+(Y-.5).^2;
        % 8层山
        for i=1:8
            % 每次都生成一次二维随机噪声，并取其中一行的数据
            h=abs(ifftn(exp(5i*rand(800))./r.^1.05)).*(CLX.*CLY).*10;
            nh=(8-i)*30+h(400,:);
            if i==1,nh=nh.*.8;end
            hm=ceil(max(nh));
            CMesh=zeros([hm,800,3]);
            
            % 颜色矩阵构造，
            tcol=colFurther+(colCloser-colFurther)./8.*(i);
            tcol=hsv2rgb(tcol./[360,100,100]);
            CMesh(:,:,1)=tcol(1);
            CMesh(:,:,2)=tcol(2);
            CMesh(:,:,3)=tcol(3);
            
            % 用nan数值框出山的轮廓
            alp=ones(hm,800);
            alp((1:hm)'>nh)=nan;
            
            % 绘制山峰
            image([-50,850],[0,hm],CMesh,'AlphaData',alp.*0.98,'Tag','SL');        
        end
    end
% =========================================================================
% 一个线性插值的渐变图生成函数
    function colorMat=vColorMat(matSize,colorList)
        yList=((0:(matSize(2)-1))./(matSize(2)-1))';
        xList=ones(1,matSize(1));
        % 线性插值
        colorMat(:,:,1)=(colorList(1,1)+yList.*(colorList(2,1)-colorList(1,1)))*xList;
        colorMat(:,:,2)=(colorList(1,2)+yList.*(colorList(2,2)-colorList(1,2)))*xList;
        colorMat(:,:,3)=(colorList(1,3)+yList.*(colorList(2,3)-colorList(1,3)))*xList;
    end
end