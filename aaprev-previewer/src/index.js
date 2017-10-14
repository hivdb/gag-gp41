import React from 'react'
import ReactDOM from 'react-dom'
import './index.css'
import Previewer from './previewer'
import registerServiceWorker from './registerServiceWorker'

ReactDOM.render(<Previewer />, document.getElementById('root'))
registerServiceWorker()
