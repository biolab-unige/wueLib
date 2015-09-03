function addCrossair(figureHandle)
%ADDCROSSAIR Description
%	 = ADDCROSSAIR(IN) Long description
%
		set(figureHandle, ...
			 'WindowButtonDownFcn', @clickFcn, ...
			 'WindowButtonUpFcn', @unclickFcn, ...
				'KeyPressFcn',@saveFcn);

		% Set up cursor text
		allLines = findobj(figureHandle, 'type', 'line');
		hText = nan(1, length(allLines));
		for id = 1:length(allLines)
			 hText(id) = text(NaN, NaN, '', ...
					'Parent', get(allLines(id), 'Parent'), ...
					'BackgroundColor', 'yellow', ...
					'Color', get(allLines(id), 'Color'));
		end

		% Set up cursor lines
		allAxes = findobj(figureHandle, 'Type', 'axes');
		hCur = nan(1, length(allAxes));
		for id = 1:length(allAxes)
			 hCur(id) = line([NaN NaN], ylim(allAxes(id)), ...
					'Color', 'black', 'Parent', allAxes(id));
		end


		function clickFcn(varargin)
				% Initiate cursor if clicked anywhere but the figure
				disp('Calling');
				if strcmpi(get(gco, 'type'), 'figure')
					 set(hCur, 'XData', [NaN NaN]);                % <-- EDIT
					 set(hText, 'Position', [NaN NaN]);            % <-- EDIT
				else
					 set(figureHandle, 'WindowButtonMotionFcn', @dragFcn)
					 pt = dragFcn();
					 set(figureHandle,'UserData',[ get(figureHandle,'UserData') pt(1) ]);
					
				end
		end

		function pt = dragFcn(varargin)
				% Get mouse location
				pt = get(gca, 'CurrentPoint');

				% Update cursor line position
				set(hCur, 'XData', [pt(1), pt(1)]);

				% Update cursor text
				for idx = 1:length(allLines)
					 xdata = get(allLines(idx), 'XData');
					 ydata = get(allLines(idx), 'YData');
					 if pt(1) >= xdata(1) && pt(1) <= xdata(end)
							y = interp1(xdata, ydata, pt(1));
							set(hText(idx), 'Position', [pt(1), y], ...
								 'String', sprintf('(%0.2f, %0.2f)', pt(1), y));
					 else
							set(hText(idx), 'Position', [NaN NaN]);
					 end
				end

		end

		function unclickFcn(varargin)
			 set(figureHandle, 'WindowButtonMotionFcn', '');
		end

		function saveFcn(varargin)
		%SAVEFCN Description
		%	 = SAVEFCN(VARARGIN) Long description
		%
		disp('Saving points');
			cursorValue = get(figureHandle,'UserData');
			assignin('base','cursorValue',cursorValue);
		end


end


