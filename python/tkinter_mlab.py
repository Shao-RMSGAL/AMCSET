import tkinter as tk
from mayavi import mlab
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, SceneEditor
from traits.api import HasTraits, Instance, on_trait_change
from traitsui.api import View, Item
from pyface.qt import QtGui, QtCore
from tvtk.pyface.api import TvtkScene
import numpy as np

import os
os.environ['ETS_TOOLKIT'] = 'qt4'
os.environ['QT_API'] = 'pyqt'

# The QWidget containing the visualization
class MayaviQWidget(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        layout.setSpacing(0)
        self.visualization = Visualization()

        # The edit_traits call will generate the widget to embed.
        self.ui = self.visualization.edit_traits(parent=self, 
                                                 kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)

class Visualization(HasTraits):
    scene = Instance(MlabSceneModel, ())

    @on_trait_change('scene.activated')
    def update_plot(self):
        x, y = np.mgrid[-2:2:100j, -2:2:100j]
        z = x * np.exp(-x**2 - y**2)
        self.scene.mlab.surf(x, y, z)

    traits_view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                            height=250, width=300, show_label=False),
                       resizable=True)

def main():
    # Don't create a new QApplication, it would unhook the Events
    # set by Traits on the existing QApplication. Simply use .instance()
    app = QtGui.QApplication.instance()
    container = QtGui.QWidget()
    container.setWindowTitle("Embedding Mayavi in a PyQt4 Application")
    # define a "complex" layout to test the behaviour
    layout = QtGui.QGridLayout(container)

    # put some stuff around mayavi
    label_list = [QtGui.QLabel(container) for i in range(3)]
    for i, label in enumerate(label_list):
        label.setText("label %d" % i)
        layout.addWidget(label, i, 0)
    mayavi_widget = MayaviQWidget(container)
    layout.addWidget(mayavi_widget, 0, 1, len(label_list), 1)
    container.show()
    window = QtGui.QMainWindow()
    window.setCentralWidget(container)
    window.show()

    # Start the main event loop.
    app.exec_()

if __name__ == "__main__":
    main()