// CustomizableCellRenderer - Modified TableCellRenderer for customized table cells

// Programmed by Yair M. Altman: altmany(at)gmail.com
// $Revision: 1.6 $  $Date: 2010/09/28 18:02:40 $

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;

public class CustomizableCellRenderer extends DefaultTableCellRenderer implements TableCellRenderer
{
	private boolean _debug = false;
	private boolean _disabled = false;
	private boolean _rowStriping = true;

	private Color _bgcolor = getBackground();
	private Color _fgcolor = getForeground();

	private final Color _disabledColor = new Color(0.925F,  0.914F,  0.847F);  // gray
	private final Color _stripeColor   = new Color(0.9608F,0.9608F, 0.9608F);  // light blue

	private final Hashtable _colorNamesHashtable = new Hashtable();

	private final Hashtable _cellBgColorHashtable   = new Hashtable();
	private final Hashtable _cellFgColorHashtable   = new Hashtable();
	private final Hashtable _cellTooltipHashtable   = new Hashtable();
	private final Hashtable _cellAlignmentHashtable = new Hashtable();
	private final Hashtable _cellFontHashtable      = new Hashtable();

	@Override
	public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column)
	{
		JComponent cell = (JComponent) super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
		if (_debug) {
			System.out.println(row + "," + column + " => value: " + value);
		}

		// Field cells may not be empty - indicate with the chosen bgcolor
		boolean emptyCell = false;
		String valueStr = "" + value;
		if ((value == null) || (valueStr.trim().length() == 0))
		{
			emptyCell = true;
			cell.setBackground(_bgcolor);
		} else if (_rowStriping) {
			Color bgcolor = (row%2 == 0) ? Color.white : _stripeColor;
			cell.setBackground(bgcolor);
		} else {
			//valueStr = (String) value;
			cell.setBackground(Color.white);
		}

		// Try to fix sorted-row reordering
		try
		{
			TableSorter model = (TableSorter) table.getModel();
			row = model.modelIndex(row);
		}
		catch (Exception ex)
		{
			// never mind...
		}

		if (_disabled)
		{
			cell.setEnabled(false);
			cell.setBackground(_bgcolor);
		}

		Vector rowColVector = getRowColVector(row, column);

		// Check whether a specific cell alignment was requested
		Integer alignment = (Integer) _cellAlignmentHashtable.get(rowColVector);
		if (alignment != null)
		{
			try
			{
				((JLabel)cell).setHorizontalAlignment(alignment.intValue());
			} catch (Exception ex) {
				alignment = null;
			}
		}
		if (alignment == null)  // also do this for caught exceptions
		{
			// Default alignment: if this is a string, align left - otherwise right
			try
			{
				java.lang.String text = (java.lang.String) value;
				if ((text.charAt(0) == '-') || (text.charAt(0) == '.') ||
						((text.charAt(0) >= '0') && (text.charAt(0) <= '9'))) {
					((JLabel)cell).setHorizontalAlignment(SwingConstants.RIGHT);
				} else {
					((JLabel)cell).setHorizontalAlignment(SwingConstants.LEFT);
				}
			}
			catch (Exception ex)
			{
				((JLabel)cell).setHorizontalAlignment(SwingConstants.LEFT);
			}
		}

		// If this cell should have a specific color, then use it
		Color cellBgColor = (Color) _cellBgColorHashtable.get(rowColVector);
		Color cellFgColor = (Color) _cellFgColorHashtable.get(rowColVector);
		if (cellBgColor != null) {
			cell.setBackground(cellBgColor);
		}
		if (cellFgColor != null) {
			cell.setForeground(cellFgColor);
		} else {
			cell.setForeground(_fgcolor);
		}

		// Check whether a specific cell font was requested
		Font font = (Font) _cellFontHashtable.get(rowColVector);
		if (font != null)
		{
			try
			{
				cell.setFont(font);
			} catch (Exception ex) {
				// never mind...
			}
		}

		// If this cell is selected, then highlight with a light-blue color
		if (isSelected) // && !(cell.getBackground().equals(Color.white)))
		{
			//cell.setBackground(cell.getBackground().darker());
			float[] rgb = cell.getBackground().getRGBComponents(null);
			rgb[0] *= .8;
			rgb[1] *= .8;  // darken the R&G components only, to highlight the blue component
			cell.setBackground(new Color(rgb[0],rgb[1],rgb[2]));
			cell.setForeground(Color.black);
		}

		// If this cell should have a specific tooltip, then use it
		String cellTooltip = (String) _cellTooltipHashtable.get(rowColVector);
		if ((cellTooltip == null) || (cellTooltip.length() == 0))
		{
			// No specific tooltip set, so use the cell's string value as the tooltip
			if (value == null) {
				cell.setToolTipText(null);
			} else if (valueStr.length() > 200)
			{
				// Split long tooltip text into several smaller lines
				String tipText = "<html>";
				int MAX_CHARS_PER_LINE = 150;
				int strLen = valueStr.length();
				for (int lineIdx=0; lineIdx <= strLen/MAX_CHARS_PER_LINE; lineIdx++) {
					tipText = tipText.concat(valueStr.substring(lineIdx*MAX_CHARS_PER_LINE,Math.min((lineIdx+1)*MAX_CHARS_PER_LINE,strLen))).concat("<br>");
				}
				cell.setToolTipText(tipText);
			} else {
				cell.setToolTipText(valueStr);
			}
		}
		else
		{
			cell.setToolTipText(cellTooltip);
		}

		return cell;
	}

	// Constructors (+ main bgcolor input arg)
	public CustomizableCellRenderer()
	{
		this(new Color(0.925F, 0.914F, 0.847F));  // gray
	}

	public CustomizableCellRenderer(float[] rgb)
	{
		this(new Color(rgb[0], rgb[1], rgb[2]));
	}

	public CustomizableCellRenderer(float r, float g, float b)
	{
		this(new Color(r,g,b));
	}

	public CustomizableCellRenderer(String colorName)
	{
		this();
		setTableBgColor(colorNameToColor(colorName));
	}

	public CustomizableCellRenderer(Color bgcolor)
	{
		super();
		//setOpaque(false);
		prepareColorNamesHash();
		setTableBgColor(bgcolor);
	}

	// Main BG color
	public void setTableBgColor(float[] rgb)
	{
		setTableBgColor(new Color(rgb[0], rgb[1], rgb[2]));
	}

	public void setTableBgColor(float r, float g, float b)
	{
		setTableBgColor(new Color(r,g,b));
	}

	public void setTableBgColor(String colorName)
	{
		setTableBgColor(colorNameToColor(colorName));
	}

	public void setTableBgColor(Color color)
	{
		_bgcolor = color;
		if (_bgcolor == null) {
			_bgcolor = Color.white;
		}
		setBackground(_bgcolor);
	}

	public Color getTableBgColor()
	{
		return _bgcolor;
	}

	// Main FG color
	public void setTableFgColor(float[] rgb)
	{
		setTableFgColor(new Color(rgb[0], rgb[1], rgb[2]));
	}

	public void setTableFgColor(float r, float g, float b)
	{
		setTableFgColor(new Color(r,g,b));
	}

	public void setTableFgColor(String colorName)
	{
		setTableFgColor(colorNameToColor(colorName));
	}

	public void setTableFgColor(Color color)
	{
		_fgcolor = color;
		if (_fgcolor == null) {
			_fgcolor = Color.black;
		}
		setForeground(_fgcolor);
	}

	public Color getTableFgColor()
	{
		return _fgcolor;
	}

	// Cell BG color
	public void resetCellBgColors()
	{
		_cellBgColorHashtable.clear();
	}

	public void setCellBgColor(int row, int column, float[] rgb)
	{
		Color color = new Color(rgb[0], rgb[1], rgb[2]);
		setCellBgColor(row, column, color);
	}

	public void setCellBgColor(int row, int column, float r, float g, float b)
	{
		Color color = new Color(r,g,b);
		setCellBgColor(row, column, color);
	}

	public void setCellBgColor(int row, int column, String colorName)
	{
		if (colorName == null)
		{
			Vector rowColVector = getRowColVector(row, column);
			_cellBgColorHashtable.remove(rowColVector);  // a null value indicates removal
		}
		else
		{
			Color color = colorNameToColor(colorName);
			setCellBgColor(row, column, color);
		}
	}

	public void setCellBgColor(int row, int column, Color color)
	{
		Vector rowColVector = getRowColVector(row, column);
		/*
		if (color == null)
			color = Color.white;                         // Hashtables cannot accept nulls...
		_cellBgColorHashtable.put(rowColVector, color);
		 */
		if (color == null) {
			_cellBgColorHashtable.remove(rowColVector);  // a null value indicates removal
		} else {
			_cellBgColorHashtable.put(rowColVector, color);
		}
	}

	public Color getCellBgColor(int row, int column)
	{
		Vector rowColVector = getRowColVector(row, column);
		return (Color) _cellBgColorHashtable.get(rowColVector);
	}

	// Cell FG color
	public void resetCellFgColors()
	{
		_cellFgColorHashtable.clear();
	}

	public void setCellFgColor(int row, int column, float[] rgb)
	{
		Color color = new Color(rgb[0], rgb[1], rgb[2]);
		setCellFgColor(row, column, color);
	}

	public void setCellFgColor(int row, int column, float r, float g, float b)
	{
		Color color = new Color(r,g,b);
		setCellFgColor(row, column, color);
	}

	public void setCellFgColor(int row, int column, String colorName)
	{
		Color color = colorNameToColor(colorName);
		setCellFgColor(row, column, color);
	}

	public void setCellFgColor(int row, int column, Color color)
	{
		Vector rowColVector = getRowColVector(row, column);
		/*
		if (color == null)
			color = Color.white;                         // Hashtables cannot accept nulls...
		_cellFgColorHashtable.put(rowColVector, color);
		 */
		if (color == null) {
			_cellFgColorHashtable.remove(rowColVector);  // a null value indicates removal
		} else {
			_cellFgColorHashtable.put(rowColVector, color);
		}
	}

	public Color getCellFgColor(int row, int column)
	{
		Vector rowColVector = getRowColVector(row, column);
		return (Color) _cellFgColorHashtable.get(rowColVector);
	}

	// Cell tooltip
	public void resetCellTooltips()
	{
		_cellTooltipHashtable.clear();
	}

	public void setCellTooltip(int row, int column, String text)
	{
		Vector rowColVector = getRowColVector(row, column);
		//System.out.println(row + "," + column + " => value: " + text);
		/*
		if (text == null)
			text = "";                                  // Hashtables cannot accept nulls...
		_cellTooltipHashtable.put(rowColVector, text);
		 */
		if (text == null) {
			_cellTooltipHashtable.remove(rowColVector);  // a null value indicates removal
		} else {
			_cellTooltipHashtable.put(rowColVector, text);
		}
	}

	public String getCellTooltip(int row, int column)
	{
		Vector rowColVector = getRowColVector(row, column);
		return (String) _cellTooltipHashtable.get(rowColVector);
	}

	// Cell alignment
	public void resetCellAlignments()
	{
		_cellAlignmentHashtable.clear();
	}

	public void setCellAlignment(int row, int column, int alignment)
	{
		Vector rowColVector = getRowColVector(row, column);
		_cellAlignmentHashtable.put(rowColVector, new Integer(alignment));
	}

	public int getCellAlignment(int row, int column)
	{
		Vector rowColVector = getRowColVector(row, column);
		Object alignment = _cellAlignmentHashtable.get(rowColVector);
		return (alignment == null) ? -1 : ((Integer)alignment).intValue();
	}

	// Cell Font
	public void resetCellFonts()
	{
		_cellFontHashtable.clear();
	}

	public void setCellFont(int row, int column, Font font)
	{
		Vector rowColVector = getRowColVector(row, column);
		if (font == null) {
			_cellFontHashtable.remove(rowColVector);  // a null value indicates removal
		} else {
			_cellFontHashtable.put(rowColVector, font);
		}
	}

	public Font getCellFont(int row, int column)
	{
		Vector rowColVector = getRowColVector(row, column);
		Font font = (Font) _cellFontHashtable.get(rowColVector);
		return (font == null) ? getFont() : font;
	}

	// Utility functions
	public void setDebug(boolean flag)
	{
		_debug = flag;
	}

	public void setDisabled(boolean flag)
	{
		_disabled = flag;
	}

	public void setRowStriping(boolean flag)
	{
		_rowStriping = flag;
	}

	private Vector getRowColVector(int row, int column)
	{
		Vector rowColVector = new Vector();
		rowColVector.addElement(new Integer(row));
		rowColVector.addElement(new Integer(column));
		return rowColVector;
	}

	private void prepareColorNamesHash()
	{
		_colorNamesHashtable.put("r",       Color.red);
		_colorNamesHashtable.put("red",     Color.red);
		_colorNamesHashtable.put("g",       Color.green);
		_colorNamesHashtable.put("green",   Color.green);
		_colorNamesHashtable.put("b",       Color.blue);
		_colorNamesHashtable.put("blue",    Color.blue);
		_colorNamesHashtable.put("c",       Color.cyan);
		_colorNamesHashtable.put("cyan",    Color.cyan);
		_colorNamesHashtable.put("m",       Color.magenta);
		_colorNamesHashtable.put("magenta", Color.magenta);
		_colorNamesHashtable.put("y",       Color.yellow);
		_colorNamesHashtable.put("yellow",  Color.yellow);
		_colorNamesHashtable.put("k",       Color.black);
		_colorNamesHashtable.put("black",   Color.black);
		_colorNamesHashtable.put("w",       Color.white);
		_colorNamesHashtable.put("white",   Color.white);
	}

	private Color colorNameToColor(String colorName)
	{
		Color color = (Color) _colorNamesHashtable.get(colorName.toLowerCase());
		return color != null ? color : _disabledColor;
	}
}
