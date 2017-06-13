function shiftAxes(handle, direction, amount)
% Shift all axes in a direction
% Directions: 'left' , 'right' , 'up' , 'down'
% shiftAxes(handle, direction, amount);
% shiftAxes(direction, amount);
if nargin<3, handle=gcf; end

switch direction

	case 'left'
	ind=1;	amount = -1 * amount;
	
	case 'right'
	ind=1;
	
	case 'up'
	ind=2;

	case 'down'
	ind=2; amount = -1 * amount;

	otherwise
	error('Direction must be : ''left'' , ''right'' , ''up'' , ''down''.');
end

children=get(handle,'Children');

	for n=1:length(children)
	if strcmpi(get(children(n),'type'),'axes')
	pos=get(children(n),'Position');
	pos(ind)=pos(ind)+amount;
	set(children(n),'Position',pos);
	end
	end

end
