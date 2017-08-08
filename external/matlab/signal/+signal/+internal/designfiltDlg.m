function choice = designfiltDlg(dlgstrs,dlgcontrolstrs)
%designfiltDlg Design filter assistant popup dialog
% choice = designfiltDlg(dlgstrs) creates a dialog assistant box with
% strings and uicontrol text provided to the function.
  
  fileNameMode = (length(dlgstrs(:)) == 5);

  dlg = signal.internal.dontshowassistantdlg;
  dlg.render;
  dlg.setDialogStrings(dlgstrs,fileNameMode)    
  dlg.setUicontrolStrings(dlgcontrolstrs) 
  
  set(dlg.Figure,'Visible','on');
  
  waitfor(dlg.Figure)
  
  choice = dlg.DlgChoice;
  
end