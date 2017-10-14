import React from 'react'
import PrevalenceDataLoader from './data-loader'
import PrevalenceDataViewer from './data-viewer'

export default class Previewer extends React.Component {
  constructor () {
    super(...arguments)
    this.state = {
      prevalenceData: null,
      wildType: null,
      regionData: null,
      prevalenceDataLoaded: false
    }
  }

  onPrevalenceDataChange (prevalenceData) {
    const prevalenceDataLoaded = prevalenceData !== null
    this.setState({prevalenceData, prevalenceDataLoaded})
  }

  onWTChange (wildType) {
    this.setState({wildType})
  }

  render () {
    const {prevalenceDataLoaded} = this.state

    return <div>
      <PrevalenceDataLoader
       onWTChange={this.onWTChange.bind(this)}
       onChange={this.onPrevalenceDataChange.bind(this)} />
      {prevalenceDataLoaded ? <PrevalenceDataViewer {...this.state} /> : null}
    </div>
  }
}
