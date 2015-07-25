% Trey Huffine, chuffine, 000274923
% EF 230 Section 1
% Project - Create an interactive golf game
% Use a figure and commands from the "player" to track distance, strength, and direction
% of a golf shot and display the path. Game should display outputs and
% strokes from the shot.
% Use user inputs (of various forms), functions, subfunctions, and Matlab
% outputs to code the game

function game = golf_trey_huffine
% Initilizes our game of golf.  Just start, click, and play!
clear; close; clc;

%% Hole 1
c=golf230('getcourse'); % get the course data from the provided file and save it to a varible
golf230('showhole',c(1)); % show the first hole
pt1 = [10 10]; % initial spot of ball. tee
[pt1, pt2] = golf230('rbline',pt1); % plot the line to get the speed of the first shot. pt1 is the initial point and pt2 is the point clicked
x1 = pt1(1); % initial x-point
y1 = pt1(2); % initial y-point
x2 = pt2(1); % final x-point of line
y2 = pt2(2); % final y point of line
speed = mag(pt1,pt2); % calculate the speed using the input line from a subfunction
angle = ang(pt1,pt2); % calculate angle using a subfunction
strokes1 = 1; % counting the first stroke

[xy, t, stat] = golf230('getpath',c(1),x1,y1,speed,angle) % saving golf230 outputs from path of shot. xy holds the xy points of the shot. t is the time of the shot. stat displays if the shot was made

% calulates these using the subfunctions below
start_distance = sd(xy); % variable for the distance from the hole before was taken
hole_distance = hd(xy); % variable to show the hole distance at the end of the shot
net_distance = nd(pt1,xy); % variable to show the net distance of the shot traveled from where the shot was taken
travel_total = tt(xy,t); % variable to show the distance traveled from where the shot was taken


% if loop to show the path of the ball from pause statements for the shots.
% the pause shows the path of the ball in increments of the shot by pausing
for i = 1:length(t);
    plot(xy(i,1),xy(i,2),'wo');
if i<length(t);
    pause(t(i+1)-t(i)); % using the pause command as a function of the time shows a rough simulation of the speed and slowing down of the ball
else
    pause(t(end)-t(end-1));
end
end


if stat == 0; % message dispalay if shot is missed ( values from above)
stats = sprintf('Strokes: %.0f \nStart distance from hole: %.2f inches \nSpeed: %.2f in/sec \nAngle: %.2f degrees \nTotal distance traveled: %.2f inches \nNet distance traveled: %.2f inches \nEnd distance from hole: %.2f inches',strokes1,start_distance,speed,angle,travel_total,net_distance,hole_distance);

game = msgbox(stats);

waitfor(game);
end


if stat == 2; % display error message if something is invalid
    wrong = warndlg('Warning Invalid Shot');
    waitfor(wrong);
end

if stat == 1; % display message box if the shot is made (values from above)
    stats = sprintf('Strokes: %.0f \nStart distance from hole: %.2f inches \nSpeed: %.2f in/sec \nAngle: %.2f degrees \nTotal distance traveled: %.2f inches \nNet distance traveled: %.2f inches \n\nYou finished the hole!',strokes1,start_distance,speed,angle,travel_total,net_distance,hole_distance);

    game = msgbox(stats);

    waitfor(game);
end

while stat == 0; % while loop run until the shot is made
    clf; % clear the last shot
    pt1 = [xy(end,1) xy(end,2)]; % set new starting point for current shot
    golf230('showhole',c(1)); % recall the course
    [pt1, pt2] = golf230('rbline',pt1); % plot the new shot. pt1 is the initial point and pt2 is the point clicked
    x1 = pt1(1); % initial x point of current shot
    y1 = pt1(2); % initial y point of current shot
    x2 = pt2(1); % final x of click
    y2 = pt2(2); % final y of click
    speed = mag(pt1,pt2); % speed of current shot
    angle = ang(pt1,pt2); % angle of current shot
    strokes1 = strokes1 + 1; % counting strokes of the shots

    [xy, t, stat] = golf230('getpath',c(1),x1,y1,speed,angle); % saving golf230 outputs to variables obtained from shot inputs. xy holds the xy points of the shot. t is the time of the shot. stat displays if the shot was made

    % calulates these using the subfunctions below
    start_distance = sd(xy); % start distance from the hole of the current shot
    hole_distance = hd(xy); % distance from the hole after the shot path has completed
    net_distance = nd(pt1,xy); % net distance that the ball travels from the starting point of the current shot
    travel_total = tt(xy,t); % total travel distance of the current shot

for i = 1:length(t); % if loop to plot shot
    plot(xy(i,1),xy(i,2),'wo');
if i<length(t);
    pause(t(i+1)-t(i));
else
    pause(t(end)-t(end-1));
end
end

if stat == 0; % message to display if shot is not made
    stats = sprintf('Strokes: %.0f \nStart distance from hole: %.2f inches \nSpeed: %.2f in/sec \nAngle: %.2f degrees \nTotal distance traveled: %.2f inches \nNet distance traveled: %.2f inches \nEnd distance from hole: %.2f inches',strokes1,start_distance,speed,angle,travel_total,net_distance,hole_distance);

    game = msgbox(stats);

    waitfor(game);
end

if stat == 1; % message to display if shot is made
    stats = sprintf('Strokes: %.0f \nStart distance from hole: %.2f inches \nSpeed: %.2f in/sec \nAngle: %.2f degrees \nTotal distance traveled: %.2f inches \nNet distance traveled: %.2f inches \n\nYou finished the hole!',strokes1,start_distance,speed,angle,travel_total,net_distance);

    game = msgbox(stats);

    waitfor(game);
end

if stat == 2; % message to display if there is an error in the shot
    warndlg('Warning Invalid Shot');
end


end

%% Hole 2
clf; % clear the first hole

golf230('showhole',c(2)); % pull of the information for hole 2
pt1 = [10 10]; % initial point of shot. tee
[pt1, pt2] = golf230('rbline',pt1); % click to get the second point to use for velocity and angle. pt1 is the initial point and pt2 is the point clicked
x1 = pt1(1); % x of initial point
y1 = pt1(2); % y of initial point
x2 = pt2(1); % x-value of click
y2 = pt2(2); % y-value of click
speed = mag(pt1,pt2); % speed calculated from click using a subprogram
angle = ang(pt1,pt2); % angle calculate from click using a sub program
strokes2 = 1; % counting the first stroke

[xy, t, stat] = golf230('getpath',c(2),x1,y1,speed,angle); % saving golf230 outputs to variables obtained from shot inputs. xy holds the xy points of the shot. t is the time of the shot. stat displays if the shot was made

% calulates these using the subfunctions below
start_distance = sd(xy); % start distance of the current shot from the hole
hole_distance = hd(xy); % distance from the hole after the shot is complete
net_distance = nd(pt1,xy); % net distance from the starting point of the current shot
travel_total = tt(xy,t); % total travel distance of the current shot

for i = 1:length(t); % if loop to see the shot path of the ball
    plot(xy(i,1),xy(i,2),'wo');
if i<length(t);
    pause(t(i+1)-t(i));
else
    pause(t(end)-t(end-1));
end
end

if stat == 0; % message to display if first shot is not made
stats = sprintf('Strokes: %.0f \nStart distance from hole: %.2f inches \nSpeed: %.2f in/sec \nAngle: %.2f degrees \nTotal distance traveled: %.2f inches \nNet distance traveled: %.2f inches \nEnd distance from hole: %.2f inches',strokes2,start_distance,speed,angle,travel_total,net_distance,hole_distance);

game = msgbox(stats);

waitfor(game);
end

if stat == 2; % error message to display if program messes up
    warndlg('Warning Invalid Shot');
end

if stat == 1; % message to display if shot is made
    stats = sprintf('Strokes: %.0f \nStart distance from hole: %.2f inches \nSpeed: %.2f in/sec \nAngle: %.2f degrees \nTotal distance traveled: %.2f inches \nNet distance traveled: %.2f inches \n\nYou finished the hole!',strokes2,start_distance,speed,angle,travel_total,net_distance);

    game = msgbox(stats);

    waitfor(game);
end

while stat == 0; % while loop to run until shot is made
    clf; % clear current shot
    pt1 = [xy(end,1) xy(end,2)]; % save the final spot of ball as new starting place.
    golf230('showhole',c(2)); % pull up the new shot
    [pt1, pt2] = golf230('rbline',pt1); % click on screen to get velocity and angle of shot. pt1 is the initial point and pt2 is the point clicked
    x1 = pt1(1); % x of initial location
    y1 = pt1(2); % y of initial location
    x2 = pt2(1); % x of point clicked
    y2 = pt2(2); % y of point clicked
    speed = mag(pt1,pt2); % subprogram used to decide speed of ball from location clicked
    angle = ang(pt1,pt2); % subprogram used to calculate angle of shot
    strokes2 = strokes2 + 1; % keep track of current strokes on current hole

    [xy, t, stat] = golf230('getpath',c(2),x1,y1,speed,angle); % saving golf230 outputs to variables obtained from shot inputs. xy holds the xy points of the shot. t is the time of the shot. stat displays if the shot was amde


    start_distance = sd(xy); % calulates these using the subfunctions below
    hole_distance = hd(xy); % all are same as stated above
    net_distance = nd(pt1,xy);
    travel_total = tt(xy,t);

    for i = 1:length(t); % if loop to plot the path of the ball as it goes through the shot
        plot(xy(i,1),xy(i,2),'wo');
        if i<length(t);
            pause(t(i+1)-t(i));
        else
            pause(t(end)-t(end-1));
        end
    end

    if stat == 0; % message to display if the shot is missed
        stats = sprintf('Strokes: %.0f \nStart distance from hole: %.2f inches \nSpeed: %.2f in/sec \nAngle: %.2f degrees \nTotal distance traveled: %.2f inches \nNet distance traveled: %.2f inches \nEnd distance from hole: %.2f inches',strokes2,start_distance,speed,angle,travel_total,net_distance,hole_distance);

        game = msgbox(stats);

        waitfor(game);
    end

    if stat == 1; % message to display if the shot is made
        stats = sprintf('Strokes: %.0f \nStart distance from hole: %.2f inches \nSpeed: %.2f in/sec \nAngle: %.2f degrees \nTotal distance traveled: %.2f inches \nNet distance traveled: %.2f inches \n\nYou finished the hole!',strokes2,start_distance,speed,angle,travel_total,net_distance);

        game = msgbox(stats);

        waitfor(game);
    end

    if stat == 2; % message to display if there is an error
        warndlg('Warning Invalid Shot');
    end

end

strokes_total = strokes1 + strokes2; % calculate total strokes from both holes

% message to display after you have finished the game and how many strokes
% it took
final = sprintf('Congratulations, you finished the course in %.0f strokes! \n\n You''re the Tiger Woods of MatLab Golf!',strokes_total);

the_end = msgbox(final);

waitfor(the_end);

close all; % close the game

return

function d = mag(pt1, pt2);
% Function to calculate mag of shot line to find speed
% is uses the initial and final point of the rbline

x1 = pt1(1); % takes the value of our initial and final points clicked and extracts the x and y values of them
x2 = pt2(1);
y1 = pt1(2);
y2 = pt2(2);

d = sqrt((x2-x1).^2+(y2-y1).^2);

return

function a = ang(pt1,pt2);
% function to calculate angle of shot using the shot line
% it uses the initial and final points of the rbline

x1 = pt1(1); % takes the value of our initial and final points clicked and extracts the x and y values of them
x2 = pt2(1);
y1 = pt1(2);
y2 = pt2(2);

a = atan2((y2-y1),(x2-x1)).*(180./pi);

return

function s = sd(xy);
% function to calculate the start distance from the hole. takes the
% magnitude of the vector from the inital place of the shot to the location of the hole

x1 = xy(1,1); % takes the initial x-value and y-value of our shot
y1 = xy(1,2);

s = sqrt((x1-200).^2+(y1-150).^2); % calculates the distance from hole using the magnitude of the distance of the hole to the initial x and y points

return

function h = hd(xy);
% function to calculate the distance from the end of the shot to the hole.
% takes the magnitude of the vector from the final spot of the ball to the
% hole

x2 = xy(end,1); % extracts the final x and y points of the shot
y2 = xy(end,2);

h = sqrt((x2-200).^2+((y2)-150).^2); % find the magnitude of the vector from the the final x and y value to the hole

return

function n = nd(pt1,xy);
% funtion to calculate the net distance traveled. finds the magnitude of the vector from
% where the shot was taken to where the ball stops

x1 = pt1(1); % extracts the initial x and y value of the shot
y1 = pt1(2);
x2 = xy(end,1); % extracts the final x and y value of the shot
y2 = xy(end,2);

n = sqrt((x2-x1).^2+(y2-y1).^2); % find the magnitude of the vector connecting the final value of the x and y of the shot to the initial value of the x and y of the shot

return

function T = tt(xy,t);
% function to calculate the total distance traveled. adds the magnitude of each consecutive
% vector of the points plotted

travel = 0 ;
for ii = 2:length(t) % finds the vector length (magnitude) for each consecutive shot and adds them together
    xx = sum(abs(xy(ii,1)-xy(ii-1,1))); % x value of the shot
    yy = sum(abs(xy(ii,2)-xy(ii-1,2))); % y value of teh shot
    travel = travel + sqrt(xx.^2+yy.^2); % magnitude of the vector using the x and y values
end

T = travel;

return
