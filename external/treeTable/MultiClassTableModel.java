// Copyright: Yair Altman
// $Revision: 1.0 $  $Date: 2013/06/20 23:25:37 $

//import java.util.Vector;

import java.util.Vector;

import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableModel;

/**
 * MultiClassTableModel is a decorator for TableModels; adding the ability
 * to set any of the JTable's columns' class to a specific Class.
 * MultiClassTableModel uses by default the value of the cell in the top
 * row to determine the column class. The user can override this using
 * a custom setColumnClass() method.
 * <p/>
 * 
 * @author Yair Altman
 * @version 1.0 06/20/13
 */

public class MultiClassTableModel extends DefaultTableModel {
	private static final long serialVersionUID = 1L;
	protected TableModel tableModel;
	private JTableHeader tableHeader;

	public MultiClassTableModel() {
		// init hashtable here
	}

	public MultiClassTableModel(TableModel tableModel) {
		this();
		setTableModel(tableModel);
	}

	public MultiClassTableModel(java.lang.Object[][] data, java.lang.Object[] headers) {
		super(data, headers);
	}

	@SuppressWarnings("rawtypes")
	public MultiClassTableModel(Vector data, Vector headers) {
		super(data, headers);
	}

	public MultiClassTableModel(TableModel tableModel, JTableHeader tableHeader) {
		this();
		setTableHeader(tableHeader);
		setTableModel(tableModel);
	}

	public TableModel getTableModel() {
		return tableModel;
	}

	public void setTableModel(TableModel tableModel) {
		this.tableModel = tableModel;
		fireTableStructureChanged();
	}

	public JTableHeader getTableHeader() {
		return tableHeader;
	}

	public void setTableHeader(JTableHeader tableHeader) {
		this.tableHeader = tableHeader;
	}

	// TableModel interface methods

	@Override
	public int getRowCount() {
		return (tableModel == null) ? super.getRowCount() : tableModel.getRowCount();
	}

	@Override
	public int getColumnCount() {
		return (tableModel == null) ? super.getColumnCount() : tableModel.getColumnCount();
	}

	@Override
	public String getColumnName(int column) {
		return (tableModel == null) ? super.getColumnName(column) : tableModel.getColumnName(column);  // fix for null case
	}

	@Override
	public Object getValueAt(int row, int column) {
		return (tableModel == null) ? super.getValueAt(row,column) : tableModel.getValueAt(row, column);  // fix for null case
	}

	/*
	// Yair 11/1/2007: several additional required interfaces
	@SuppressWarnings("rawtypes")
	@Override
	public void addRow(Vector vec) {
		if (tableModel != null) {
			tableModel.addRow(vec);
		} else {
			super.addRow(vec);
		}
	}
	@Override
	public void addRow(Object[] vec) {
		try {
			if (tableModel != null) {
				tableModel.addRow(vec);
			} else {
				super.addRow(vec);
			}
		} catch (Exception e) {
			// do nothing
		}
	}

	@Override
	public void addColumn(Object name) {
		if (tableModel != null) {
			tableModel.addColumn(name);
		} else {
			super.addColumn(name);
		}
	}
	@SuppressWarnings("rawtypes")
	@Override
	public void addColumn(Object name, Vector data) {
		if (tableModel != null) {
			tableModel.addColumn(name,data);
		} else {
			super.addColumn(name,data);
		}
	}
	@Override
	public void addColumn(Object name, Object[] data) {
		if (tableModel != null) {
			tableModel.addColumn(name,data);
		} else {
			super.addColumn(name,data);
		}
	}

	@SuppressWarnings("rawtypes")
	@Override
	public void insertRow(int beforeRow, Vector data) {
		if (tableModel != null) {
			tableModel.insertRow(beforeRow, data);
		} else {
			super.insertRow(beforeRow, data);
		}
	}
	@Override
	public void insertRow(int beforeRow, Object[] data) {
		if (tableModel != null) {
			tableModel.insertRow(beforeRow, data);
		} else {
			super.insertRow(beforeRow, data);
		}
	}

	@Override
	public void removeRow(int rowIdx) {
		if (tableModel != null) {
			tableModel.removeRow(rowIdx);
		} else {
			super.removeRow(rowIdx);
		}
	}

	@Override
	public void moveRow(int a, int b, int c) {
		if (tableModel != null) {
			tableModel.moveRow(a,b,c);
		} else {
			super.moveRow(a,b,c);
		}
	}

	@Override
	public void setNumRows(int numRows) {
		if (tableModel != null) {
			tableModel.setNumRows(numRows);
		} else {
			super.setNumRows(numRows);
		}
	}
	@Override
	public void setRowCount(int numRows) {
		try {
			if (tableModel != null) {
				tableModel.setRowCount(numRows);
			} else {
				super.setRowCount(numRows);
			}
		} catch (Exception e) {
			// do nothing
		}
	}
	@Override
	public void setColumnCount(int numCols) {
		if (tableModel != null) {
			tableModel.setColumnCount(numCols);
		} else {
			super.setColumnCount(numCols);
		}
	}

	@SuppressWarnings("rawtypes")
	@Override
	public void setColumnIdentifiers(Vector colNames) {
		if (tableModel != null) {
			tableModel.setColumnIdentifiers(colNames);
		} else {
			super.setColumnIdentifiers(colNames);
		}
	}
	@Override
	public void setColumnIdentifiers(Object[] colNames) {
		try {
			if (tableModel != null) {
				tableModel.setColumnIdentifiers(colNames);
			} else {
				super.setColumnIdentifiers(colNames);
			}
		} catch (Exception e) {
			// do nothing
		}
	}

	@SuppressWarnings("rawtypes")
	@Override
	public Vector getDataVector() {
		return (tableModel == null) ? super.getDataVector() : tableModel.getDataVector();
	}
	@SuppressWarnings("rawtypes")
	@Override
	public void setDataVector(Vector data1, Vector data2) {
		if (tableModel != null) {
			tableModel.setDataVector(data1, data2);
		} else {
			super.setDataVector(data1, data2);
		}
	}
	@Override
	public void setDataVector(Object[][] data1, Object[] data2) {
		if (tableModel != null) {
			tableModel.setDataVector(data1, data2);
		} else {
			super.setDataVector(data1, data2);
		}
	}

	public String getColumnModel(int column) {
		return (tableModel == null) ? super.getColumnName(column) : tableModel.getColumnName(column);  // fix for null case
	}
	 */

	// Yair: infer the column class from the actual data.
	// Otherwise, LEXICAL_COMPARATOR (not COMPARABLE_COMPARATOR) will be used also for numbers
	@SuppressWarnings({ "unchecked", "rawtypes" })
	@Override
	public Class getColumnClass(int column) {
		Class objClass;
		try {
			//return tableModel.getColumnClass(column);    // Yair - old code
			if (getRowCount() > 0) {
				objClass = getValueAt(0, column).getClass();
			}
			else {
				objClass = (tableModel == null) ? super.getColumnClass(column) : tableModel.getColumnClass(column);  // revert back to model's class if no data; also-nullfix
			}
		}
		catch (Exception e)
		{
			objClass = "".getClass();  // null fix
		}
		//System.out.println("MultiClassTableModel.getColumnClass("+column+") = " + objClass);
		return objClass;
	}
}
