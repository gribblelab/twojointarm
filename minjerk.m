function [T,H,Hd,Hdd] = minjerk(H1,H2,t,n)

% Given hand initial position H1=(x1,y1), final position H2=(x2,y2) and movement duration t,
% and the total number of desired sampled points n,
% Calculates the hand path H over time T that satisfies minimum-jerk.
% Also returns derivatives Hd and Hdd. From the paper:
%   Flash T. & Hogan N. (1985) "The coordination of arm
%   movements: an experimentally confirmed mathematical model."
%   Journal of Neuroscience 5(7): 1688-1703.

T = linspace(0,t,n)';
H = zeros(n,2);
Hd = zeros(n,2);
Hdd = zeros(n,2);
H(:,1) = H1(1) + ((H1(1)-H2(1))*(15*((T/t).^4) - (6*(T/t).^5) - (10*(T/t).^3)));
H(:,2) = H1(2) + ((H1(2)-H2(2))*(15*((T/t).^4) - (6*(T/t).^5) - (10*(T/t).^3)));
Hd(:,1) = (H1(1) - H2(1))*(-30*(T.^4/t^5) + 60*(T.^3/t^4) - 30*(T.^2/t^3));
Hd(:,2) = (H1(2) - H2(2))*(-30*(T.^4/t^5) + 60*(T.^3/t^4) - 30*(T.^2/t^3));
Hdd(:,1) = (H1(1) - H2(1))*(-120*(T.^3/t^5) + 180*(T.^2/t^4) - 60*(T/t^3));
Hdd(:,2) = (H1(2) - H2(2))*(-120*(T.^3/t^5) + 180*(T.^2/t^4) - 60*(T/t^3));

end
