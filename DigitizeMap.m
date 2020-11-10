% ----------------------------------------------------------------------
% Digitize a map
% The map should be of Mercator Projection with North direction upward and
% uniform latitude. 
% Collected data will be stored as gmt psxy input file. 
%
% Author: Yi LUO
% Built: 2019-09-11  11:00
% Last Modified: 2019-09-17  13:30
% ----------------------------------------------------------------------

%% input parameters --------------------------------------------
% mapfile='Hall2000.PNG';
mapfile='Hall2000_2.PNG';
% mapfile='Ray1.jpg';

% Anchor points which should be far enough in two dimension
% anchor=[100,10;140,-10];
anchor=[120 2;124 -6];

% output_name=[mapfile,'.test','.xy'];
output_name='hall2000_sulawesi_suture.txt';
% output_name='ray.txt';

%% Main Part ---------------------------------------------------
% ----------------------------------------------------------------------
if exist(output_name,'file')
    warning('"%s" already exists!',output_name);
    prompt='Type 1 to add; 0 to overwrite; Ctrl+C to escape:';
    flag1=input(prompt);
    if flag1==1
        output=fopen(output_name,'a');
    elseif flag1==0
        output=fopen(output_name,'w');
    end
else
    output=fopen(output_name,'w');
end

prompt='Type 1 for curves(psxy); 2 for points(psxy,pstext); 3 for polygons(psxy):';
purpose=input(prompt);

oldanchorflag=exist('anchor_x','var');
if oldanchorflag==1
    prompt='Anchor points exist! Enter to heritage; 0 to overwrite:';
    oldanchorflag2=input(prompt);
    if ~isempty(oldanchorflag2)
        anchorflag=1;
    else
        anchorflag=0;
    end
else
    anchorflag=1;
end

y=imread(mapfile);
figure
imshow(y);
hold on

if anchorflag==1
disp('Click the Northwest and Southeast anchor points with coordinates of:');
disp(num2str(anchor));
[anchor_x,anchor_y]=ginput(2);
end

%% Curves ------------------------------------------------------------
if purpose==1
    
    flag=1;
    while flag~=0
        disp('Click to draw a curve. Press enter to end the curve. Press again to finish work.')
% ------ faster but no display while drawing
%         [temp_x,temp_y]=ginput;
%         
%         if isempty(temp_x)
%             break
%         end
%         x=(temp_x-anchor_x(1))./(anchor_x(2)-anchor_x(1)).*(anchor(2,1)-anchor(1,1))+anchor(1,1);
%         y=(temp_y-anchor_y(2))./(anchor_y(1)-anchor_y(2)).*(anchor(1,2)-anchor(2,2))+anchor(2,2);
%         fprintf(output,'%f %f\n',[x';y']);
%         fprintf(output,'>\n');
%         plot(temp_x,temp_y,'-x','color','r')
% ------ slower but show every point clicked
        previous_x=zeros(0);
        previous_y=previous_x;
        flag2=1;
        while flag2~=0
            [temp_x,temp_y]=ginput(1);    
            if isempty(temp_x)
                break;
            end
            flag2=2;
            x=(temp_x-anchor_x(1))./(anchor_x(2)-anchor_x(1)).*(anchor(2,1)-anchor(1,1))+anchor(1,1);
            y=(temp_y-anchor_y(2))./(anchor_y(1)-anchor_y(2)).*(anchor(1,2)-anchor(2,2))+anchor(2,2);
            fprintf(output,'%f %f\n',[x,y]);
            disp([x,y])
            plot([previous_x,temp_x],[previous_y,temp_y],'-x','color','r')
            previous_x=temp_x;
            previous_y=temp_y;
        end
        if flag2==1
            break
        end
        fprintf(output,'>\n'); 
% ------        
        disp('Recorded')
    end
    
%% Points, texts -------------------------------------------------------
elseif purpose==2
    flag=1;
    while flag~=0
        disp('Click to draw a point. Press enter to finish work.')
        [temp_x,temp_y]=ginput(1);
        
        if isempty(temp_x)
            break
        end
        x=(temp_x-anchor_x(1))./(anchor_x(2)-anchor_x(1)).*(anchor(2,1)-anchor(1,1))+anchor(1,1);
        y=(temp_y-anchor_y(2))./(anchor_y(1)-anchor_y(2)).*(anchor(1,2)-anchor(2,2))+anchor(2,2);
        fprintf(output,'%f %f',[x,y]);
        plot(temp_x,temp_y,'o','color','r')
        disp([x,y])
        prompt='Type the text here, or leave it blank:';
        param_str=input(prompt,'s');
        fprintf(output,' %s\n',param_str);
    end
    
%% Polygons ------------------------------------------------------
elseif purpose==3
        flag=1;
    while flag~=0
        disp('Click to draw a polygon. Press enter to end the polygon. Press again to finish work.')
        previous_x=zeros(0);
        previous_y=previous_x;
        flag2=1;
        n=1;
        while flag2~=0
            [temp_x,temp_y]=ginput(1);    
            if isempty(temp_x)
                break;
            end
            flag2=2;
            x=(temp_x-anchor_x(1))./(anchor_x(2)-anchor_x(1)).*(anchor(2,1)-anchor(1,1))+anchor(1,1);
            y=(temp_y-anchor_y(2))./(anchor_y(1)-anchor_y(2)).*(anchor(1,2)-anchor(2,2))+anchor(2,2);
            fprintf(output,'%f %f\n',[x,y]);
            disp([x,y])
            plot([previous_x,temp_x],[previous_y,temp_y],'-x','color','r')
            previous_x=temp_x;
            previous_y=temp_y;
            if n==1
                tempx1=temp_x;
                tempy1=temp_y;
                x1=x;
                y1=y;
                n=2;
            end
        end
        if flag2==1
            break
        end
        fprintf(output,'%f %f\n',[x1,y1]);
        disp([x1,y1])
        plot([previous_x,tempx1],[previous_y,tempy1],'-x','color','r')
        
        fprintf(output,'>\n'); 
% ------        
        disp('Recorded')
    end
    %% end -------------------------------------------------------
end
fclose(output);
disp('Finished!')

