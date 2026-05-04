#!/usr/bin/env python3
import sys
from PySide6 import QtWidgets
from src.astrotouch.gui import AstroTouchWindow

def main():
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("Fusion")
    
    window = AstroTouchWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
